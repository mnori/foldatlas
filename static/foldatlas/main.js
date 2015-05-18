
// Define the website's controller. All page changes must go through this badboy
function BrowserController(config) {

	var instance = this;
	
	this.nucsPerRow = 80;
	this.staticBase = config.staticBaseUrl;
	this.reactivities = {}

	// initialise: function(config) {
	// 	this.staticBase = window.browserControllerConfig;
	// },

	this.init = function() {
		instance.drawReactivities();
	}

	// Jump to a specific transcript page
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

	// HTML5 change URL method
	this.changeUrl = function(title, url) {
	    if (typeof (history.pushState) != "undefined") {
	        var obj = { Page: title, Url: url };
	        history.pushState(obj, obj.Page, obj.Url);
	    } else {
	        alert("Your browser does not support HTML5. Please upgrade it.");
	    }
	}

	// Visualises the reactivity data.
	this.drawReactivities = function(reactivities) {

		var data = this.getReactivitiesJson()
		if (data == null) {
			return;
		}

		console.log(data);

		var nDataRows = data.length;
		var nChartRows = Math.ceil(nDataRows / this.nucsPerRow);

		console.log("nChartRows: "+nChartRows);

		// Define chart dimensions including axis panelMargins
		var panelMargin = {top: 15, right: 60, bottom: 30, left: 70}
		var panelTotDims = {x: 898, y: 100}

		// dims without margins
		var panelDims = { 
			x: panelTotDims.x - panelMargin.left - panelMargin.right,
			y: panelTotDims.y - panelMargin.bottom - panelMargin.top
		}

		// Total dimensions of chart across all panels and margins.
		var chartDims = {
			x: panelTotDims.x,
			y: panelTotDims.y * nChartRows
		}

		// var	width = panelTotDims.x - panelMargin.left - panelMargin.right,
		// 	height = panelTotDims.y - panelMargin.bottom - panelMargin.top;

		// Init the chart's container element
		var chart = d3.select("#reactivities-chart")
			.attr("width", chartDims.x)
			.attr("height", chartDims.y)

		// Define the scales
		var yScale = d3.scale.linear()

			// range maps to the pixel dimensions.

		    // domain describes the range of data to show.
		    .range([panelDims.y, 0])
		    .domain([0, d3.max(data, function(d) { return d.reactivity; })]);

		// when there is no reactivity data, degrade gracefully
		if (isNaN(yScale.domain()[1])) { 
			
			chart.append("text")
		      .attr("transform", "translate("+(panelTotDims.x / 2)+", "+(panelTotDims.y / 2)+")")
		      .style("text-anchor", "middle")
		      .text("No reactivity data");
			return;
		}
		
		for (var rowN = 0; rowN < nChartRows; rowN++) { // each iteration = 1 chart row
			var start = rowN * this.nucsPerRow;
			var end = start + this.nucsPerRow;

			if (end > nDataRows) {
				end = nDataRows;
			}

			var dataSlice = data.slice(start, end);

			var nucsThisRow = end - start;

			// this is for panel positioning.
			panelYOffset = rowN * panelTotDims.y;

			// Shows nucleotide numbers
			var rangeX = parseInt(panelDims.x * ((nucsThisRow - 1) / this.nucsPerRow));

			var xScale = d3.scale.linear()
			    .range([0, rangeX], .1) 
			    .domain([start - 0.5, (end - 1) + 0.5])

		   	// Create axis objects
			var xAxis = d3.svg.axis()
			    .scale(xScale)
			    .orient("bottom")
			    .ticks(nucsThisRow)
				.tickFormat(function(d, i) { return data[d].nuc; })

			// xAxis.select("tick")
			// 	.append("rect")
			// 	.setAttr("width", 10)
			// 	.setAttr("height", 10);

			var yAxis = d3.svg.axis()
			    .scale(yScale)
			    .orient("left")
	    		.ticks(3); // how many ticks to show.

			// Add x-axis objects to the chart. X axis ticks have a background,
			// which depends on the data value.
			var bgWidth = parseInt(panelDims.x / this.nucsPerRow) + 1;

			chart.append("g")
				.attr("class", "x axis")
				.attr("transform", "translate("+
					panelMargin.left+","+
					(panelYOffset + panelDims.y + panelMargin.top)
				+")")
				.call(xAxis)
				.selectAll(".tick")
				.insert("rect", ":first-child")
				.attr("transform", "translate("+(-bgWidth / 2)+", "+10+")")
				.attr("width", bgWidth)
				.attr("height", 10)

				// highlight nucleotides with missing reactivities
				.attr("class", function(n, i) {
					return (dataSlice[i].reactivity == null) ? 
						"missing-bg" : "not-missing-bg";
				})

			// Add y-axis objects to the chart
			chart.append("g")
				.attr("class", "y axis")
				.attr("transform", "translate("+
					panelMargin.left+","+
					(panelYOffset + panelMargin.top)+
				")")
				.call(yAxis);

			// Add x axis label
		    chart.append("text")
		        .attr("transform", "translate("+
		        	(panelMargin.left + (panelDims.x / 2))+","+
		        	(panelYOffset + panelDims.y + panelMargin.top)+
	        	")")
		        .style("text-anchor", "middle")
		        .attr("dy", "2.7em");
		        // .text("Nucleotide");

			// Add y-axis label
			chart.append("text")
		        .attr("transform", "rotate(-90)")
		        .attr("y", panelMargin.left) // this is actually X direction, because we rotated.
		        .attr("x", (-panelYOffset - (panelMargin.top + (panelDims.y / 2))))
		        .attr("dy", "-2.7em")
		        .style("text-anchor", "middle")
		        .text("Reactivity");

	        // Add length label
	        var panelDimsX = panelDims.x * (nucsThisRow / this.nucsPerRow);
			chart.append("text")
		        .attr("x", panelMargin.left + panelDimsX)
		        .attr("y", panelYOffset + panelMargin.top + panelDims.y)
		        .style("text-anchor", "left")
		        .attr("dy", "1.3em")
		        .text(end);

		    // add the actual line chart
			var lineGen = d3.svg.line()
			    .x(function(d) { return panelMargin.left + xScale(d.position); })
			    .y(function(d) { return panelYOffset + panelMargin.top + yScale(d.reactivity); });

	 		chart.append("path")
				.datum(dataSlice) // get data specific to this row
				.attr("class", "line")
				.attr("d", lineGen);

			// var barWidth = panelDims.x / this.nucsPerRow;

			// bar = chart.selectAll("g").select("g").data(data).enter()
			// 	.attr("transform", function(d, i) { 
			// 		console.log(i);
			// 		return "translate(" + (
			// 			panelMargin.left + (i * barWidth)) + ", "+
			// 			panelMargin.top+
			// 		")"; 
				// });

			// Add rectangles to the bars
			// // // Add bar elements to the chart
			// var bar = chart.selectAll("g")
			// 	.data(data.slice(start, end))
			// 	.enter().append("g")
			// 	.attr("transform", function(d, i) { 
			// 		console.log("pos: "+d.position);
			// 		// console.log("i"+ i);
			// 		return "translate(" + (panelMargin.left + (i * barWidth)) + ", "+panelMargin.top+")"; 
			// 	})

			// bar.append("rect")
			// 	.attr("y", function(d) { return yScale(d.reactivity); })
			// 	.attr("height", function(d) { return panelDims.y - yScale(d.reactivity); })
			// 	.attr("width", barWidth);	

		} // end looping through chart rows

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
		// 	.attr("transform", function(d, i) { return "translate(" + (panelMargin.left + (i * barWidth)) + ", "+panelMargin.top+")"; })

		// // Add rectangles to the bars
		// bar.append("rect")
		// 	.attr("y", function(d) { return yScale(d.reactivity); })
		// 	.attr("height", function(d) { return height - yScale(d.reactivity); })
		// 	.attr("width", barWidth);	
	}

	this.getReactivitiesJson = function() {
		var html = $("#reactivities-json").html();

		if (html == undefined) { // this means no reactivity data to show
			return null;
		}
		var json = $.parseJSON(html);
		return json;
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