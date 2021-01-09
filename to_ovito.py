import ovito 

# want to import the files into ovito 
pipeline = ovito.io.import_file("Snapshots/snapshot_t*.xyz", 
    columns = ["Position.X", "Position.Y", "Position.Z", "Particle Type"])
# Now we want to change the radius of cells in the simulation
# Get the particle data to manipulate it:
data = pipeline.source.data
particle_info = data.particles.particle_types.types
# Now change the radius - start with setting it to 2 for all cells
# Obviously radii are different in vivo, but for illustration we try at 2
for ptype in particle_info:
    ptype.radius = 2.0 
# TODO figure out why some cells are different sizes?
# Now load the data into the renderer
pipeline.add_to_scene()

vp = ovito.vis.Viewport(camera_dir = (1, 1, -1)) # find that ortho view is too jittery
vp.zoom_all()

# set the background colour of the animation here, using RGB format
background_colour = [0.239, 0.239, 0.239] # default grey is 0.239, 0.239, 0.239 (higher for lighter)
# Set the resolution here - put it to 1080p for default
resolution = [1920, 1080] # this is the standard 1920x1080 pixels (x,y)
# Now render the animation and save it as you wish - can be .mp4, .mov, .avi, .gif 
file_name = "animation.avi"
# Now render the animation - can adjust the fps here if you wish
# TODO choose a best renderer - ovito gives us two options
vp.render_anim(file_name, size=resolution, fps=20, background=background_colour)

# Works!