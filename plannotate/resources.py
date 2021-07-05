import pkg_resources

def get_resource(group, name):
    return pkg_resources.resource_filename(__package__, f"data/{group}/{name}")

def get_image(name):
    return get_resource("images", name)

def get_template(name):
    return get_resource("templates", name)

def get_example_fastas():
    return get_resource("fastas", "")
