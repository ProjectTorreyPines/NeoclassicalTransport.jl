using Documenter, NeoclassicalTransport

# Call functions
open(joinpath(@__DIR__, "src/api.md"), "w") do f
    println(f, "# API Reference\n")
    for page in keys(NeoclassicalTransport.document)
        if page == :Expressions
            continue
        end
        println(f, "## $page\n")
        println(f, "```@docs")
        for item in NeoclassicalTransport.document[page]
            println(f, "$item")
        end
        println(f, "```")
    end
end

makedocs(;
    modules=[NeoclassicalTransport],
    format=Documenter.HTML(;analytics="G-65D8V8C8VQ"),
    sitename="NeoclassicalTransport",
    checkdocs=:none,
    pages=["index.md", "api.md", "License" => "license.md", "Notice" => "notice.md"],
    warnonly=true
)

# Deploy docs
# This function deploys the documentation to the gh-pages branch of the repository.
# The main documentation that will be hosted on
# https://projecttorreypines.github.io/NEO.jl/stable
# will be built from latest release tagged with a version number.
# The development documentation that will be hosted on
# https://projecttorreypines.github.io/NEO.jl/dev
# will be built from the latest commit on the chosen devbranch argument below.
# For testing purposes, the devbranch argument can be set to WIP branch like "docs".
# While merging with master, the devbranch argument should be set to "master".
deploydocs(;
    repo="github.com/ProjectTorreyPines/NeoclassicalTransport.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "v#.#"]
)
