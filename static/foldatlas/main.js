
// Define the website's controller. All page changes must go through this badboy
function BrowserController(config) {

	this.staticBase = config.staticBaseUrl;

	// initialise: function(config) {
	// 	this.staticBase = window.browserControllerConfig;
	// },

	this.selectTranscript = function(transcriptID) {
		this.changeUrl(transcriptID, "/transcript/"+transcriptID)
		$("#transcript-data").html("Loading...")
		$("#transcript-data").load("/ajax/transcript/"+transcriptID)
	}

	this.changeUrl = function(title, url) {
	    if (typeof (history.pushState) != "undefined") {
	        var obj = { Page: title, Url: url };
	        history.pushState(obj, obj.Page, obj.Url);
	    } else {
	        alert("Your browser does not support HTML5. Please upgrade it.");
	    }
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