#!/usr/bin/env bash
# update.sh — unattended refresh of the package data from the rule bank; meant for cron.
#
#   build/update.sh              rebuild; if a rule changed, test, commit and push
#   build/update.sh --install    add the hourly cron entry for this script (idempotent)
#   build/update.sh --remove     take the cron entry out again
#
# What one run does:
#   1. skips if the working tree has uncommitted changes (somebody is working here);
#   2. git pull --ff-only, then build/build_data.jl, which replaces data/ only when the new
#      data are no worse (exit 3: regression, exit 4: credit change — both leave data/ alone
#      and raise a desktop notification, because they need a person);
#   3. if files changed: runs the test suite; on success commits and pushes, on failure
#      restores the tree and notifies;
#   4. carries the data (and RULES.md, FORMAT.md, NOTICE.md, the README's credits block) into the
#      sibling Python and R packages (../quadriceps-py,
#      ../quadriceps-r) when they exist: copy, run their tests, commit, push.
#
# Run it on ONE machine only (it needs the project's sync folder; the bank is the same
# everywhere, and two machines would race to push the same commit).
# Log: ~/.local/state/quadriceps/update.log

set -u
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1   # stay out of the way
PKG=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SYNC=${QUADRICEPS_SYNC:-$HOME/Dropbox/oldDesignedQuadrature-sync}
STATE=$HOME/.local/state/quadriceps; mkdir -p "$STATE"
LOG=$STATE/update.log
CRONLINE="43 * * * * $PKG/build/update.sh >/dev/null 2>&1 # Quadriceps.jl data update"

case "${1:-}" in
  --install) { crontab -l 2>/dev/null | grep -vF "# Quadriceps.jl data update"; echo "$CRONLINE"; } | crontab -
             echo "installed: $CRONLINE"; exit 0 ;;
  --remove)  crontab -l 2>/dev/null | grep -vF "# Quadriceps.jl data update" | crontab -; echo "removed"; exit 0 ;;
esac

exec 9>"$STATE/lock"; flock -n 9 || exit 0
exec >>"$LOG" 2>&1
say()    { echo "$(date '+%F %T') $*"; }
notify() { say "NOTIFY: $*"
           DBUS_SESSION_BUS_ADDRESS=${DBUS_SESSION_BUS_ADDRESS:-unix:path=/run/user/$(id -u)/bus} \
             notify-send -u critical "Quadriceps update" "$*" 2>/dev/null || true; }

# Released Julia only: the project's launcher picks the newest 1.13.x; fall back to juliaup's default.
J=$(bash "$SYNC/project/julia/symq_julia.sh" 2>/dev/null | tail -1)
[ -x "${J:-}" ] || J=$HOME/.juliaup/bin/julia
[ -x "$J" ] || { notify "no julia found"; exit 1; }
[ -d "$SYNC/project/julia/rules" ] || { say "no rule bank at $SYNC; skipped"; exit 0; }

clean() { [ -z "$(git -C "$1" status --porcelain)" ]; }

cd "$PKG" || exit 1
clean "$PKG" || { say "working tree not clean; skipped"; exit 0; }
git pull -q --ff-only || { notify "git pull failed in $PKG"; exit 1; }

OUT=$(nice -n 19 "$J" --project=build build/build_data.jl "$SYNC" 2>&1); RC=$?
case $RC in
  0) ;;
  3|4) notify "$(echo "$OUT" | grep -A3 -E '^(REGRESSION|CREDIT CHANGE)' | head -4)"; git checkout -q -- . ; exit $RC ;;
  *) notify "build_data.jl failed (exit $RC): $(echo "$OUT" | tail -2)"; git checkout -q -- . ; git clean -qfd data; exit $RC ;;
esac

if ! clean "$PKG"; then
  CHANGES=$(echo "$OUT" | sed -n 's/^CHANGE: //p')
  if nice -n 19 "$J" --project=. -e 'using Pkg; Pkg.test()' >"$STATE/test.log" 2>&1; then
    git add -A
    git commit -q -m "Data update from the rule bank" -m "${CHANGES:-generated files refreshed}"
    git push -q && say "pushed: $(echo "$CHANGES" | tr '\n' ';')" || notify "git push failed in $PKG"
  else
    notify "tests failed after a data rebuild; tree restored (see $STATE/test.log)"
    git checkout -q -- . ; git clean -qfd data; exit 1
  fi
fi

# Sibling packages: same data, their own tests.  sibling <dir> <data dir inside it> <test command…>
sibling() {
  local dir=$1 data=$2; shift 2
  [ -d "$dir/.git" ] || return 0
  clean "$dir" || { say "$dir: working tree not clean; skipped"; return 0; }
  git -C "$dir" pull -q --ff-only || { notify "git pull failed in $dir"; return 1; }
  mkdir -p "$dir/$data"
  rsync -a --delete --exclude rules128.bin --exclude index128.tsv "$PKG/data/" "$dir/$data/"   # quadruple precision: Julia only
  sed 's|(credits.md)|(NOTICE.md)|' "$PKG/docs/src/rules.md" > "$dir/RULES.md"
  sed -e 's|(credits.md)|(NOTICE.md)|' -e 's|^# Data format|# Data format of rules.bin|' "$PKG/docs/src/format.md" > "$dir/FORMAT.md"
  cp "$PKG/NOTICE.md" "$dir/NOTICE.md"
  # the GENERATED credits block of the twin's README, from this package's README
  awk -v src="$PKG/README.md" '
    BEGIN { while ((getline l < src) > 0) { if (l ~ /BEGIN GENERATED credits/) { on = 1; continue }
                                            if (l ~ /END GENERATED credits/) on = 0
                                            if (on) block = block l "\n" } }
    /BEGIN GENERATED credits/ { print; printf "%s", block; skip = 1; next }
    /END GENERATED credits/   { skip = 0 }
    !skip' "$dir/README.md" > "$dir/README.md.new" && mv "$dir/README.md.new" "$dir/README.md"
  clean "$dir" && return 0
  if ( cd "$dir" && nice -n 19 "$@" ) >"$STATE/test-$(basename "$dir").log" 2>&1; then
    git -C "$dir" add -A
    git -C "$dir" commit -q -m "Data update from Quadriceps.jl" -m "$(git -C "$PKG" log -1 --format=%b)"
    git -C "$dir" push -q && say "$dir: pushed" || notify "git push failed in $dir"
  else
    notify "$(basename "$dir"): tests failed after a data update; tree restored"
    git -C "$dir" checkout -q -- . ; git -C "$dir" clean -qfd
  fi
}
sibling "$PKG/../quadriceps-py" src/quadriceps/data python3 -m pytest -q
sibling "$PKG/../quadriceps-r"  inst/extdata        Rscript tests/run_tests.R
say "done"
