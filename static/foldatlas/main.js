
// Define the website's controller. All page changes must go through this badboy
function BrowserController(config) {

	this.staticBase = config.staticBaseUrl;

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
    window.browserController = new BrowserController(window.foldatlasConfig)
});

$("#chromosome-selector").change(function() {
	// TODO reload the page .. difficult / impossible to change the chromosome
	// ... 
});