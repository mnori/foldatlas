
// Define the website's controller. All page changes must go through this badboy
function BrowserController(config) {

	this.staticBase = config.staticBase;

	// initialise: function(config) {
	// 	this.staticBase = window.browserControllerConfig;
	// },

	this.selectTranscript = function(transcriptID) {

		// console.log(transcriptID)
		$("#transcript-data").html("Loading...")

		var ajaxUrl = "/ajax/transcript/"+transcriptID;
		$("#transcript-data").load(ajaxUrl)

		
		// // alert("ajaxUrl: "+ajaxUrl)
		// $.ajax({
		// 	url: ajaxUrl
		// }).done(function() {
		// 	alert("It worked!")
		// })
		// .. fire off ajax request
	}
}

$(document).ready(function () { 
    window.genoverse = new Genoverse(window.genoverseConfig); 
    window.browserController = new BrowserController(window.browserControllerConfig)
});

$("#chromosome-selector").change(function() {
	// TODO reload the page

	// window.genoverse.resetConfig()
	// window.genoverse.loadGenome()
	// window.genoverse.loadPlugins()
	// window.genoverse.addTracks()
	// window.genoverse.destroy()
	// window.genoverse.init()

	delete(window.genoverse)
	window.genoverse = new Genoverse(window.genoverseConfig); 
	// window.genoverse.reload(window.genoverseConfig)
});