
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
		var margin = {top: 20, right: 20, bottom: 60, left: 70}

		var totWidth = 898,
			totHeight = 248;

		var	width = totWidth - margin.left - margin.right,
			height = totHeight - margin.bottom - margin.top;

		// Init the chart's container element
		var chart = d3.select("#reactivities-chart")
			.attr("width", width + margin.left + margin.right)
			.attr("height", height + margin.bottom + margin.top)

		// Define the scales
		var yScale = d3.scale.linear()

			// range maps to the pixel dimensions.
		    .range([height, 0])

		    // domain describes the range of data to show.
		    .domain([0, d3.max(data, function(d) { return d.reactivity; })]);

		// when there is no reactivity data, degrade gracefully
		if (isNaN(yScale.domain()[1])) { 
			
			chart.append("text")
		      .attr("transform", "translate("+(totWidth / 2)+", "+(totHeight / 2)+")")
		      .style("text-anchor", "middle")
		      .text("No reactivity data");
			return;
		}
		
		var xScale = d3.scale.linear()
		    .range([0, width], .1) 
		    .domain([0, d3.max(data, function(d) { return d.position; })]);
		    // .domain([0, 80]);

	   	// Create axis objects
		var xAxis = d3.svg.axis()
		    .scale(xScale)
		    .orient("bottom")
		    .ticks(10);

		var yAxis = d3.svg.axis()
		    .scale(yScale)
		    .orient("left")
    		.ticks(5); // how many ticks to show.

		// Add x-axis objectsto the chart
		chart.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate("+margin.left+"," + (height + margin.top) + ")")
			.call(xAxis);

		// Add x axis label
	    chart.append("text")
	        .attr("transform", "translate(" + (margin.left + (width / 2)) + " ," + (height + margin.top) + ")")
	        .style("text-anchor", "middle")
	        .attr("dy", "2.7em")
	        .text("Nucleotide");

		// Add y-axis objects to the chart
		chart.append("g")
			.attr("class", "y axis")
			.attr("transform", "translate("+margin.left+", "+margin.top+")")
			.call(yAxis)

		// Add y-axis label
		chart.append("text")
	        .attr("transform", "rotate(-90)")
	        .attr("y", margin.left) // this is actually X direction, because we rotated.
	        .attr("x", 0 - (margin.top + (height / 2)))
	        .attr("dy", "-2.7em")
	        .style("text-anchor", "middle")
	        .text("Normalised Reactivity");

	    // add the actual line chart
		var lineGen = d3.svg.line()
		    .x(function(d) { return margin.left + xScale(d.position); })
		    .y(function(d) { return margin.top + yScale(d.reactivity); });

 		chart.append("path")
			.datum(data)
			.attr("class", "line")
			.attr("d", lineGen);

		// old code for generating bar chart
		// chart.append('svg:path')
		// 	.attr('d', lineGen(data))
		// 	.attr('stroke', 'green')
		// 	.attr('stroke-width', 2)
		// 	.attr('fill', 'none');

		// // // Define bar width, depends on n values and also width of canvas
		// var barWidth = width / data.length;

		// // Add bar elements to the chart
		// var bar = chart.selectAll("g")
		// 	.data(data)
		// 	.enter().append("g")
		// 	.attr("transform", function(d, i) { return "translate(" + (margin.left + (i * barWidth)) + ", "+margin.top+")"; })

		// // Add rectangles to the bars
		// bar.append("rect")
		// 	.attr("y", function(d) { return yScale(d.reactivity); })
		// 	.attr("height", function(d) { return height - yScale(d.reactivity); })
		// 	.attr("width", barWidth);	

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