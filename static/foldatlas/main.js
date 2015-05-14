
// Define the website's controller. All page changes must go through this badboy
function BrowserController(config) {

	var instance = this;

	this.staticBase = config.staticBaseUrl;
	this.reactivities = {}

	// initialise: function(config) {
	// 	this.staticBase = window.browserControllerConfig;
	// },

	this.selectTranscript = function(transcriptID) {
		$("#transcript-loading").show()
		this.changeUrl(transcriptID, "/transcript/"+transcriptID)
		$("#transcript-data").load("/ajax/transcript/"+transcriptID, function() {
			$("#transcript-loading").hide();

			// not the best way to do it...
			var json = $.parseJSON($("#reactivities-json").html())

			console.log(json)

			if (json.length != 0) {
				instance.drawReactivities(json);
			} else {
				alert("No reactivity data found");
			}
		})
	}

	this.changeUrl = function(title, url) {
	    if (typeof (history.pushState) != "undefined") {
	        var obj = { Page: title, Url: url };
	        history.pushState(obj, obj.Page, obj.Url);
	    } else {
	        alert("Your browser does not support HTML5. Please upgrade it.");
	    }
	}

	this.drawReactivities = function(reactivities) {

		var data = reactivities

		// Define chart dimensions including axis margins
		var margin = {left: 80, bottom: 80}
		var width = 900 - margin.left, height = 250 - margin.bottom;

		// Define the scales
		var x = d3.scale.linear()
		    .range([0, width]);

		var y = d3.scale.linear()
		    .range([height, 0]);

	   	// Create axis objects
		var xAxis = d3.svg.axis()
		    .scale(x)
		    .orient("bottom")
		    .ticks(10);

		var yAxis = d3.svg.axis()
		    .scale(y)
		    .orient("left")
    		.ticks(10);

		// Define bar width, depends on n values and also width of canvas
		var barWidth = width / data.length;
		console.log("barWidth: "+barWidth);

		// Get the chart's container element
		var chart = d3.select("#reactivities-chart")
			.attr("width", width + margin.left)
			.attr("height", height + margin.bottom)

		// Add x-axis objects to the chart
		chart.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate("+margin.left+"," + height + ")")
			.call(xAxis);

		// Add y-axis objects to the chart
		chart.append("g")
			.attr("class", "y axis")
			.attr("transform", "translate("+margin.left+", 0)")
			.call(yAxis)

		// Add bar elements to the chart
		var bar = chart.selectAll("g")
			.data(data)
			.enter().append("g")
			.attr("transform", function(d, i) { return "translate(" + (margin.left + (i * barWidth)) + ",0)"; })

		// Add rectangles to the bars
		bar.append("rect")
			.attr("y", function(d) { return y(d.reactivity); })
			.attr("height", function(d) { return height - y(d.reactivity); })
			.attr("width", barWidth);	

		// DONT add text to the bars
		// bar.append("text")
		// 	.attr("x", barWidth / 2)
		// 	.attr("y", function(d) { return y(d.reactivity) + 3; })
		// 	.attr("dy", ".75em")
		// 	.text(function(d) { return d.reactivity; });
		
		// function type(d) {
		//   d.reactivity = +d.reactivity; // coerce to number
		//   return d;
		// }
	}
}

$(document).ready(function () { 
    window.genoverse = new Genoverse(window.genoverseConfig); 
    window.browserController = new BrowserController(window.foldatlasConfig);
});

$("#chromosome-selector").change(function() {
	// TODO reload the page .. difficult / impossible to change the chromosome
	// ... 
});