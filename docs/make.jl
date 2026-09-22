using Documenter, Quadriceps

makedocs(
    sitename = "Quadriceps.jl",
    modules = [Quadriceps],
    authors = "Joris Pinkse",
    format = Documenter.HTML(prettyurls = false, size_threshold = nothing),
    pages = [
        "Home" => "index.md",
        "Guide" => "guide.md",
        "Stored rules" => "rules.md",
        "Data format" => "format.md",
        "Reference" => "api.md",
        "Credits" => "credits.md",
    ],
    checkdocs = :exports,
)

deploydocs(repo = "github.com/NittanyLion/Quadriceps.jl.git", devbranch = "main", push_preview = false)
