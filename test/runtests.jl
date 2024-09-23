using Test
using Aqua

using ReefGuideAPI

@testset "Aqua" begin
    Aqua.test_unbound_args(ReefGuideAPI)
    Aqua.test_undefined_exports(ReefGuideAPI)
    Aqua.test_project_extras(ReefGuideAPI)
    Aqua.test_stale_deps(ReefGuideAPI)
end
