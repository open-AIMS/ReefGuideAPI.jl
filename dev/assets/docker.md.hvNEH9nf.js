import{_ as i,c as s,o as a,ai as t}from"./chunks/framework.3Uj6wpql.js";const c=JSON.parse('{"title":"Docker build and run","description":"","frontmatter":{},"headers":[],"relativePath":"docker.md","filePath":"docker.md","lastUpdated":null}'),l={name:"docker.md"};function n(o,e,r,h,d,u){return a(),s("div",null,e[0]||(e[0]=[t(`<h1 id="Docker-build-and-run" tabindex="-1">Docker build and run <a class="header-anchor" href="#Docker-build-and-run" aria-label="Permalink to &quot;Docker build and run {#Docker-build-and-run}&quot;">​</a></h1><p>The ReefGuideAPI.jl package has an associated <code>Dockerfile</code> and build/publish process. This means you can run an instance of the ReefGuideAPI.jl package without needing to compile/build it with a local <code>Julia</code> installation. You will be able to view the latest published versions of the Docker image on the repository packages page.</p><h2 id="Dockerfile-Configuration" tabindex="-1">Dockerfile Configuration <a class="header-anchor" href="#Dockerfile-Configuration" aria-label="Permalink to &quot;Dockerfile Configuration {#Dockerfile-Configuration}&quot;">​</a></h2><p>The Julia version is specified by the <code>JULIA_VERSION</code> arg. The full version is specified to maintain build stability, but should be bumped to the latest version of Julia when a release is published.</p><h2 id="Publish-Release" tabindex="-1">Publish Release <a class="header-anchor" href="#Publish-Release" aria-label="Permalink to &quot;Publish Release {#Publish-Release}&quot;">​</a></h2><ol><li><p>Bump version in Project.toml and <a href="../../.github/workflows/PublishDockerImage.yml">PublishDockerImage.yml</a></p></li><li><p>Create PR, merge to <code>main</code> branch</p></li><li><p>Publish Release on GitHub (this triggers the PublishDockerImage workflow)</p></li></ol><h2 id="A-note-about-MKL_jll" tabindex="-1">A note about MKL_jll <a class="header-anchor" href="#A-note-about-MKL_jll" aria-label="Permalink to &quot;A note about MKL_jll {#A-note-about-MKL_jll}&quot;">​</a></h2><p>Due to how Julia (particularly v1.11) handles precompilation, it significantly reduces the build time by explicitly installed MKL_jll before installing any of explicit project dependencies.</p><p>For this reason, the Dockerfile extracts the MKL_jll version from the Manifest file using Pkg.dependency(), precompiles this in an anonymous project, then compiles the main dependencies. This cuts the build time from around 15 minutes down to around 6-7.</p><h2 id="Mounting-files-and-required-data" tabindex="-1">Mounting files and required data <a class="header-anchor" href="#Mounting-files-and-required-data" aria-label="Permalink to &quot;Mounting files and required data {#Mounting-files-and-required-data}&quot;">​</a></h2><p>As mentioned in <a href="./@ref getting_started">Getting Started</a>, the <code>ReefGuideAPI.jl</code> package currently requires</p><ul><li><p>a <code>config.toml</code> file and</p></li><li><p>a set of input data files</p></li></ul><p>Please include these in a folder called <code>data</code> in your working directory.</p><p>When running the below commands, it is assumed you have <code>data</code> available locally with the required files.</p><p><strong>Note</strong>: Due to how Docker excludes <code>.</code> files, we have named the config file <code>config.toml</code> in the data folder. This is required to launch the server.</p><h2 id="To-build-from-src-files-using-Docker" tabindex="-1">To build from src files using Docker <a class="header-anchor" href="#To-build-from-src-files-using-Docker" aria-label="Permalink to &quot;To build from src files using Docker {#To-build-from-src-files-using-Docker}&quot;">​</a></h2><div class="language-bash vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">bash</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">docker</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> build</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> .</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> --target</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> reefguide-src</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> -t</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> reefguide</span></span></code></pre></div><h2 id="To-build-from-src-files-using-Docker-Compose" tabindex="-1">To build from src files using Docker Compose <a class="header-anchor" href="#To-build-from-src-files-using-Docker-Compose" aria-label="Permalink to &quot;To build from src files using Docker Compose {#To-build-from-src-files-using-Docker-Compose}&quot;">​</a></h2><div class="language-bash vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">bash</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">docker</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> compose</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> up</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> --build</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> reefguide-src</span></span></code></pre></div><h2 id="To-run-with-mounted-files-(launch-server)-using-Docker" tabindex="-1">To run with mounted files (launch server) using Docker <a class="header-anchor" href="#To-run-with-mounted-files-(launch-server)-using-Docker" aria-label="Permalink to &quot;To run with mounted files (launch server) using Docker {#To-run-with-mounted-files-(launch-server)-using-Docker}&quot;">​</a></h2><div class="language-bash vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">bash</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">docker</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> run</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> -p</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> 8000:8000</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> -v</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> ./data:/data/reefguide</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> reefguide</span></span></code></pre></div><h2 id="To-run-with-mounted-files-(launch-server)-using-Docker-Compose" tabindex="-1">To run with mounted files (launch server) using Docker Compose <a class="header-anchor" href="#To-run-with-mounted-files-(launch-server)-using-Docker-Compose" aria-label="Permalink to &quot;To run with mounted files (launch server) using Docker Compose {#To-run-with-mounted-files-(launch-server)-using-Docker-Compose}&quot;">​</a></h2><div class="language-bash vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">bash</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">docker</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> compose</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> up</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> reefguide-src</span></span></code></pre></div><h2 id="To-run-with-mounted-files-(interactive-shell)-using-Docker" tabindex="-1">To run with mounted files (interactive shell) using Docker <a class="header-anchor" href="#To-run-with-mounted-files-(interactive-shell)-using-Docker" aria-label="Permalink to &quot;To run with mounted files (interactive shell) using Docker {#To-run-with-mounted-files-(interactive-shell)-using-Docker}&quot;">​</a></h2><p>This will start a Julia shell where <code>ReefGuideAPI</code> is compiled and ready for use e.g.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ReefGuideAPI</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ReefGuideAPI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">start_server</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;/data/reefguide/config.toml&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language-bash vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">bash</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">docker</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> run</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> --rm</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> --interactive</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> --entrypoint=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;julia&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> --tty</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> -v</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> ./data:/data/reefguide</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> reefguide</span></span></code></pre></div>`,27)]))}const k=i(l,[["render",n]]);export{c as __pageData,k as default};
