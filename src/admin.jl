
function setup_admin_routes(config)
    @get "/health" function ()
        return json(Dict(:status => "healthy"))
    end
end
