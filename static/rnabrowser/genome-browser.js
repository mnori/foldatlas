$(document).ready(function () { 
    window.genoverse = new Genoverse(window.genoverseConfig); 
});


$("#chromosome-selector").change(function() {
	// window.genoverse.resetConfig()
	// window.genoverse.loadGenome()
	// window.genoverse.loadPlugins()
	// window.genoverse.addTracks()
	// window.genoverse.destroy()
	// window.genoverse.init()


	window.genoverse.reload(window.genoverseConfig)
});