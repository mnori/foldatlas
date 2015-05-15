
// Define the website's controller. All page changes must go through this badboy
function BrowserController(config) {

	var instance = this;

	this.staticBase = config.staticBaseUrl;
	this.reactivities = {}

	// initialise: function(config) {
	// 	this.staticBase = window.browserControllerConfig;
	// },

	this.init = function() {
		instance.drawReactivities();
	}

	this.selectTranscript = function(transcriptID) {
		$("#loading-indicator").show()
		this.changeUrl(transcriptID, "/transcript/"+transcriptID)

		$.ajax({
			url: "/ajax/transcript/"+transcriptID,
			context: this
		}).done(function(results) {
			$("#loading-indicator").hide();
			$("#transcript-data").empty();
			$("#transcript-data").html(results);
			instance.drawReactivities();
		});
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

		var data = this.getReactivitiesJson()

		// Define chart dimensions including axis margins
		var margin = {top: 20, right: 20, bottom: 40, left: 40}

		var totWidth = 898,
			totHeight = 248;

		var	width = totWidth - margin.left - margin.right,
			height = totHeight - margin.bottom - margin.top;

		// Init the chart's container element
		var chart = d3.select("#reactivities-chart")
			.attr("width", width + margin.left + margin.right)
			.attr("height", height + margin.bottom + margin.top)

		// If there's no data, show a message. Otherwise continue plotting the data.
		// ...

		var y = d3.scale.linear()
		    .range([height, 0])
		    .domain([0, d3.max(data, function(d) { return d.reactivity; })]);

		console.log(y.domain());
		console.log(y.domain()[1]);
		if (isNaN(y.domain()[1])) { // happens when there is no reactivity data
			chart.append("text")
		      .attr("transform", "translate("+(totWidth / 2)+", "+(totHeight / 2)+")")
		      .style("text-anchor", "middle")
		      .text("No reactivity data");
		   return;
		}

		// Define the scales
		var x = d3.scale.linear()
			// range maps to the pixel dimensions.
		    .range([0, width], .1) 

		    // domain describes the range of data to show.
		    .domain([0, d3.max(data, function(d) { return d.position; })]);

		



	   	// Create axis objects
		var xAxis = d3.svg.axis()
		    .scale(x)
		    .orient("bottom")
		    .ticks(10);

		var yAxis = d3.svg.axis()
		    .scale(y)
		    .orient("left")
    		.ticks(5); // how many ticks to show.

		// Define bar width, depends on n values and also width of canvas
		var barWidth = width / data.length;

		// Add x-axis objects to the chart
		chart.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate("+margin.left+"," + (height + margin.top) + ")")
			.call(xAxis);

		// Add y-axis objects to the chart
		chart.append("g")
			.attr("class", "y axis")
			.attr("transform", "translate("+margin.left+", "+margin.top+")")
			.call(yAxis)

		// Add bar elements to the chart
		var bar = chart.selectAll("g")
			.data(data)
			.enter().append("g")
			.attr("transform", function(d, i) { return "translate(" + (margin.left + (i * barWidth)) + ", "+margin.top+")"; })

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

	this.getReactivitiesJson = function() {
		var json = $.parseJSON($("#reactivities-json").html())
		return json
	}

	this.init()
}

$(document).ready(function () { 
    window.genoverse = new Genoverse(window.genoverseConfig); 
    window.browserController = new BrowserController(window.foldatlasConfig);
});

$("#chromosome-selector").change(function() {
	// TODO reload the page .. difficult / impossible to change the chromosome 
	// ... 
});