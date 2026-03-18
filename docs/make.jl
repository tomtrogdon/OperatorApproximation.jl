using Documenter, OperatorApproximation

makedocs(
    sitename = "OperatorApproximation.jl",
    authors = "Tom Trogdon <trogdon@uw.edu>",
    repo = Remotes.GitHub("tomtrogdon", "OperatorApproximation.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://tomtrogdon.github.io/OperatorApproximation.jl",
    ),
    modules = [OperatorApproximation],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Domains" => "domains.md",
        "Bases" => "bases.md",
        "Operators" => "operators.md",
        "Solving Equations" => "solvers.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :none,
)

deploydocs(
    repo = "github.com/tomtrogdon/OperatorApproximation.jl",
    devbranch = "main",
)
