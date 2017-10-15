* Define a simulation class which knows how to read and execute user-defined run files.  Use a json format, such as:
{
	"version": "1.0",
	"name": "my simulation",
	"onset_flow": [1.0,0.0],
	"airfoils": [{"type":"naca", "properties":{"code":"0012"}}],
	"motion": blah,
	"resolution": 50,
	"num_cycles": 5,
	"max_time": 3.0,
	"save_frames_to": "output/frames",
	"save_movie": "output/movie.mp4"
}

* Define a __main__ method in the simulation class, which takes as input a path/URL to the run json definition; if not provided, read it from stdin, so that the simulation can be used in UNIX terminal pipelines.
