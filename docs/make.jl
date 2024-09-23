using
    Documenter,
    DocumenterVitepress,
    DocumenterTools
using Logging

using ReefGuideAPI

Logging.disable_logging(Logging.Warn)

makedocs(;
    sitename="ReefGuideAPI.jl",
    modules=[ReefGuideAPI],
    clean=true,
    doctest=true,
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/open-AIMS/ReefGuideAPI.jl",
        devbranch="main",
        devurl = "dev"
    ),
    draft=false,
    source="src",
    build="build",
    warnonly=true
)

# Enable logging to console again
Logging.disable_logging(Logging.BelowMinLevel)

deploydocs(;
    repo="github.com/open-AIMS/ReefGuideAPI.jl.git",
    target="build", # this is where Vitepress stores its output
    branch = "gh-pages",
    devbranch="main",
    push_preview=false
)
