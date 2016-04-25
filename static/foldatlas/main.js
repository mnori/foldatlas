/**
 * Front end code for FoldAtlas.
 * @author Matthew Norris <matthew.norris@jic.ac.uk> 2015
 */

 window.faGlobals = {
	selectedTid: null,
	selectedSid: null
}

// Create the browser controller class
var BrowserController = Class.extend({

	// Constructor
	init: function(config) {
		
		this.nucsPerRow = 80;
		this.staticBase = config.staticBaseUrl;
		this.reactivities = {};
		this.drawTranscriptData();
		this.searchController = new SearchController(this);

		$("#title").click($.proxy(function() { this.jumpTo("/"); }, this));
		$("#help-button").click($.proxy(function(ev) { 
			ev.preventDefault();
			this.jumpTo("/help"); 
		}, this));

		// Detect back/forward buttons
		// We must react by changing the page for each type of URL
		$(window).on("popstate", $.proxy(function(event) {
			var url = document.location;
			this.onBackForward(url);
		}, this));

		this.initHelpLinks();

	},

	// URLS //////////////////////////////////////////////////////////

	// This countains all the dynamic URL mappings to JS and titles
	routeUrl: function(urlIn) {
		urlIn = ""+urlIn;

		if (urlIn.indexOf("#") != -1) { // internal # link - ignore
			return;
		}

		if (urlIn.indexOf("/search") != -1) {
			this.showHideLoading(); // just for consistency
			this.searchController.show();
			document.title = "FoldAtlas: Search";

		} else if (urlIn.indexOf("/transcript/") != -1) {
			var transcriptID = this.getTranscriptID(urlIn);
			this.selectTranscript(transcriptID);
			document.title = "FoldAtlas: "+transcriptID;

		} else if (urlIn.indexOf("/help") != -1) {
			this.selectHelp();
			document.title = "FoldAtlas: Help"

		} else { // url not recognised - assume we need to go home
			this.showHideLoading(); // just for consistency
			this.goHome();
			document.title = "FoldAtlas";
		}
	},

	selectHelp: function() {
		this.showLoading();
		$.ajax({
			url: "/ajax/help",
			context: this

		}).done(function(html) {
			// TODO sort out this mess
			$("#help").show()
			$("#search").hide()
			$("#d3nome").hide();
			$("#help").html(html);
			$("#transcript-data").empty();
			this.hideLoading();
			this.initHelpLinks();
		});
	},

	initHelpLinks: function() {
		// Attach listeners to each element
		$(".ref-link").each($.proxy(function(key, element) {
			var element = $(element)
			element.click(function() {

				// Highlight reference when link is clicked.
				// Remove highlighting from all other references.
				var id = element.attr("href").replace("#", "");
				$(".ref-link").each($.proxy(function(key, element) {
					var element = $(element);
					element.removeClass("highlighted");
				}, this))
				$("#"+id).addClass("highlighted");
			})
		}, this));
	},

	getTranscriptID: function(urlIn) {
		if (typeof(urlIn) === "undefined") {
			urlIn = ""+document.location;
		}
		return urlIn.split("/").pop();
	},

	onBackForward: function(url) {
		this.routeUrl(url);
	},

	// Call this within the JS to change the page.
	jumpTo: function(url) {
		this.routeUrl(url); // execute JS for this URL
		this.changeUrl(url); // change the URL bar
	},

	// HTML5 change URL method
	changeUrl: function(url, title) {
	    if (typeof (history.pushState) != "undefined") {
	    	// the title was already changed in routeUrl so we can just fetch it out here.
	        var obj = {
	        	Page: document.title,
	        	Url: url
	        };
	        history.pushState(obj, obj.Page, obj.Url);
	    } else {
	        alert("Your browser does not support HTML5. Please upgrade it.");
	    }
	},

	// /URLS /////////////////////////////////////////////////////////

	showLoading: function() {
		$("#loading-indicator").show()
	},

	hideLoading: function() {
		$("#loading-indicator").fadeOut(300);
	},

	showHideLoading: function() {
		this.showLoading();
		this.hideLoading();
	},

	// Reset to landing page
	goHome: function() {
		// TODO sort out this mess
		$("#help").hide();
		$("#search").hide();
		$("#d3nome").show();
		$("#transcript-data").empty();
		this.hideLoading();
	},

	getJsonFromElement: function(elementID) {
		var html = $("#"+elementID).html();
		if (html == undefined) { // this means no measurement data to show
			return null;
		}
		var json = $.parseJSON(html);
		return json;
	},

	// Jump to a specific transcript page
	selectTranscript: function(transcriptID) {
		window.faGlobals.selectedTid = transcriptID;
		this.showLoading();
		// this.changeUrl("/transcript/"+transcriptID, transcriptID)
		$.ajax({
			url: "/ajax/transcript/"+transcriptID,
			context: this
		}).done(function(results) {
			$("#help").hide()
			$("#search").hide()
			$("#d3nome").show();
			$("#transcript-data").empty();
			$("#transcript-data").html(results);
			this.drawTranscriptData();
			this.hideLoading();
		});
	},

	/**
	 * Draw transcript data, if available.
	 */
	drawTranscriptData: function() {
		// 1) obtain transcript coordinates
		var transcriptData = this.getJsonFromElement("transcript-json");
		if (!transcriptData) {
			return;
		}

		// 2) use those to move the brush to the right place
		var chrID = transcriptData["chromosome_id"];
		var diff = transcriptData["end"] - transcriptData["start"];
		var start = transcriptData["start"] - diff;
		var end = transcriptData["end"] + diff;

		// figure out chrInd
		var d3nome = window.d3nomeObject;
		d3nome.jumpToPosition(d3nome.chrIDToInd(chrID), [start, end], true);

		var structureData = this.getJsonFromElement("structure-json")

		// TODO better way to detect whether there is stuctural data
		if (structureData) {
			this.structureExplorer = new StructureExplorer(this)
		}

		var measurementData = this.getJsonFromElement("nucleotide-measurement-json")
		if (measurementData) {
			var transcriptID = this.getTranscriptID();
			this.drawNucleotideMeasurements(measurementData[1], transcriptID);

			// disable ribosome profiling .. for now
			// this.drawNucleotideMeasurements(measurementData[2], transcriptID);
		}
	},

	// this should be a seperate class.
	drawNucleotideMeasurements: function(experimentData, transcriptID) {
		var detailedID = "nucleotide-measurements-overview_"+experimentData["id"];
		var detailedContainerID = detailedID+"_container"
		var overviewID = "nucleotide-measurements-detailed_"+experimentData["id"];

		var moreDetailID = "more-detail_"+experimentData["id"];
		var lessDetailID = "less-detail_"+experimentData["id"];

		// Add HTML elements (TODO use template instead?)

		if (experimentData["empty"]) {
			var buf = 
				"<h2 class=\"bar\">"+experimentData["description"]+"</h2>"+
				"<div class=\"message\">No data available.</div>";
			$("#nucleotide-measurement-charts").append(buf)

		} else {
			buf = 
				"<h2 class=\"bar\">"+
					experimentData["description"]+
					"<a href=\"/download/raw_measurements/"+experimentData["id"]+"/"+this.getTranscriptID()+"\" target=\"_blank\" class=\"button download r\">"+
						"<i class=\"fa fa-download\"></i> Download raw"+
					"</a>"+
					"<a href=\"/download/measurements/"+experimentData["id"]+"/"+this.getTranscriptID()+"\" target=\"_blank\" class=\"button download\">"+
						"<i class=\"fa fa-download\"></i> Download normalised"+
					"</a>"+
				"</h2>"+
				"<div id=\""+overviewID+"_container\" class=\"nm-container\">"+
					"<a href=\"javascript:void(0)\" id=\"nm-overview-dl-button\" class=\"button svg\">"+
						"<i class=\"fa fa-file-image-o\"></i>"+
					"</a>"+
					"<svg id=\""+overviewID+"\"></svg>"+
				"</div>"+
				"<a href=\"#\" id=\""+moreDetailID+"\" class=\"nucleotide-detail button\">"+
					"<i class=\"fa fa-arrow-circle-down\"></i>&nbsp;"+
					"More detail"+
				"</a>"+
				"<a href=\"#\" id=\""+lessDetailID+"\" class=\"nucleotide-detail button\" style=\"display: none;\">"+
				"<i class=\"fa fa-arrow-circle-up\"></i>&nbsp;"+
					"Less detail"+
				"</a>"+
				"<div id=\""+detailedContainerID+"\" style=\"display: none;\" class=\"nm-container\">"+
					"<a href=\"javascript:void(0)\" id=\"nm-detailed-dl-button\" class=\"button svg\">"+
						"<i class=\"fa fa-file-image-o\"></i>"+
					"</a>"+
					"<svg id=\""+detailedID+"\"></svg>"+
				"</div>";
			
			$("#nucleotide-measurement-charts").append(buf)

			// Add button event handlers
			$("#"+moreDetailID).click($.proxy(function(ev) {
				ev.preventDefault();
				$("#"+detailedContainerID).show();
				$("#"+moreDetailID).hide();
				$("#"+lessDetailID).show();
			}, this))

			$("#"+lessDetailID).click($.proxy(function(ev) {
				ev.preventDefault();
				$("#"+detailedContainerID).hide();
				$("#"+moreDetailID).show();
				$("#"+lessDetailID).hide();
			}, this))

			// Draw the charts
			this.drawNucleotideMeasurementsOverview(overviewID, experimentData);
			this.drawNucleotideMeasurementsDetailed(detailedID, experimentData);

			var dlID = window.faGlobals.selectedTid;

			new SvgDownloader(overviewID, "nm-overview-dl-button", "reacts_"+dlID+".svg");
			new SvgDownloader(detailedID, "nm-detailed-dl-button", "reacts-detailed_"+dlID+".svg");
		}
	},

	// Draw a graph with just 1 row, showing all of the nucleotide measurements
	// on the same row. Could potentially make this zoomable for more awesomeness
	drawNucleotideMeasurementsOverview: function(svgID, experimentData) {

		// var yLabelText = (experimentData["type"] == "dms_reactivity") 
		// 	? "Reactivity"
		// 	: "Occupancy";

		var data = experimentData["data"]
		if (data == null) { // can happen
			return;
		}

		var nDataRows = data.length;

		// Define chart dimensions including axis panelMargins
		var panelMargin = {top: 15, right: 60, bottom: 30, left: 35}
		var panelTotDims = {x: 898, y: 100}

		// for showing length string
		var lengthOffset = 10

		// dims without margins
		var panelDims = { 
			x: panelTotDims.x - panelMargin.left - panelMargin.right,
			y: panelTotDims.y - panelMargin.bottom - panelMargin.top
		}

		// Init the chart's container element
		var chart = d3.select("#"+svgID)
		var chartContainer = d3.select("#"+svgID+"_container")
		var styleStr = 
			"width: "+panelTotDims.x+"px; "+
			"height: "+panelTotDims.y+"px; ";
		chart.attr("style", styleStr);
		chartContainer.attr("style", styleStr);

		var maxY = d3.max(data, function(d) { return d.measurement; });
		var maxX = data.length;

		// Define the scales
		var yScale = d3.scale.linear()
			// range maps to the pixel dimensions.
		    // domain describes the range of data to show.
		    .range([panelDims.y, 0])
		    .domain([0, maxY]);

		var xScale = d3.scale.linear()
			    .range([0, panelDims.x], .1) // what does the .1 do? Add extra space?
			    .domain([-0.5, (maxX - 1) + 0.5])

	   	// Create axis objects
		var xAxis = d3.svg.axis()
		    .scale(xScale)
		    .orient("bottom")
		    .ticks(10)
			// .tickFormat(function(d, i) { return data[d].nuc; }) // tells it to use the nucleotide letter

		var yAxis = d3.svg.axis()
		    .scale(yScale)
		    .orient("left")
    		.ticks(3); // how many ticks to show.

    	// need to add a new x axis tick for this thing.

		// Add y-axis objects to the chart
		var yAxisEl = chart.append("g")
			.attr("class", "y axis")
			.attr("transform", "translate("+panelMargin.left+","+panelMargin.top+")")
			.call(yAxis)

		chart.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate("+
				panelMargin.left+","+
				(panelMargin.top + panelDims.y)+
			")")
			.call(xAxis);

		// Add inline style to lines and paths
		var lineStyle = {	
			"fill": "none",
			"stroke": "#000",
			"stroke-width": "1px"
		}
		chart.selectAll("path").style(lineStyle);
  		chart.selectAll("line").style(lineStyle)
  			// .style({"stroke", "#000"});

		// Add y-axis label
		// chart.append("text")
	 //        .attr("transform", "rotate(-90)")
	 //        .attr("y", panelMargin.left) // this is actually X direction, because we rotated.
	 //        .attr("x", panelMargin.top - panelDims.y)
	 //        .attr("dy", "-2.7em")
	 //        .style("text-anchor", "middle")
	 //        .text(yLabelText);

        // Add length label
        var panelDimsX = panelDims.x;
		chart.append("text")
	        .attr("x", panelMargin.left + panelDimsX + lengthOffset)
	        .attr("y", panelMargin.top + panelDims.y)
	        .style("text-anchor", "left")
	        .attr("dy", "1.3em")
	        .text(data.length);

		// Draw bars
		var barWidth = parseInt(panelDims.x / data.length);
		if (barWidth < 1) {
			barWidth = 1;
		}
		var bar = chart
			.selectAll("g.nucleotide-measurement-bar")
			.data(data).enter()
			.append("g")
			.attr("class", "nucleotide-measurement-bar")
			.attr("transform", function(d) { 
				return "translate("+
					(panelMargin.left + xScale(d.position) - (barWidth / 2))+","+
					(panelMargin.top + yScale(d.measurement))+
				")";
			});

		bar.append("rect")
			.attr("height", function(d) { 
				return yScale(maxY - d.measurement);
			})
			.attr("width", barWidth)
			.style("fill", "#c33");
	},

	// Visualises the measurement data.
	drawNucleotideMeasurementsDetailed: function(svgID, experimentData) {

		// for showing length string
		var lengthOffset = 10

		// var yLabelText = (experimentData["type"] == "dms_reactivity") 
		// 	? "Reactivity"
		// 	: "Occupancy";

		var data = experimentData["data"]
		if (data == null) { // can happen
			return;
		}

		var nDataRows = data.length;
		var nChartRows = Math.ceil(nDataRows / this.nucsPerRow);

		// Define chart dimensions including axis panelMargins
		var panelMargin = {top: 15, right: 60, bottom: 30, left: 35}
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
		var chart = d3.select("#"+svgID)
			.attr("width", chartDims.x)
			.attr("height", chartDims.y)

		var maxY = d3.max(data, function(d) { return d.measurement; });

		// Define the scales
		var yScale = d3.scale.linear()
			// range maps to the pixel dimensions.
		    // domain describes the range of data to show.
		    .range([panelDims.y, 0])
		    .domain([0, maxY]);

		// when there is no measurement data, degrade gracefully
		// TODO get rid of this - handle it higher up
		if (isNaN(yScale.domain()[1])) { 
			
			chart.append("text")
		      .attr("transform", "translate("+(panelTotDims.x / 2)+", "+(panelTotDims.y / 2)+")")
		      .style("text-anchor", "middle")
		      .text("No measurement data");
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
			var rangeX = parseInt(panelDims.x * (nucsThisRow / this.nucsPerRow));

			var xScale = d3.scale.linear()
			    .range([0, rangeX], .1) 
			    .domain([start - 0.5, (end - 1) + 0.5])

		   	// Create axis objects
			var xAxis = d3.svg.axis()
			    .scale(xScale)
			    .orient("bottom")
			    .ticks(nucsThisRow)
				.tickFormat(function(d, i) { return data[d].nuc; }) // tells it to use the nucleotide letter

			var yAxis = d3.svg.axis()
			    .scale(yScale)
			    .orient("left")
	    		.ticks(3); // how many ticks to show.

			// Add x-axis objects to the chart.
			var bgWidth = parseInt(panelDims.x / this.nucsPerRow) + 1;

			var xAxisElement = chart.append("g")
				.attr("class", "x axis")
				.attr("transform", "translate("+
					panelMargin.left+","+
					(panelYOffset + panelDims.y + panelMargin.top)
				+")")
				.call(xAxis)

			xAxisElement
				.selectAll(".tick") // select the X axis tick element

				// aadd a rect to it, first child means it draws in the background
				.insert("rect", ":first-child")

				// set position and dimensions of rectangle element.
				.attr("transform", "translate("+(-bgWidth / 2)+", "+10+")")
				.attr("width", bgWidth)
				.attr("height", 10)
				.style("fill", function(n, i) { 
					var nucMissing = dataSlice[i].measurement == null
					return nucMissing ? "#bbb" : "#fff"; 
				})

			xAxisElement.selectAll("text")
				.style("font-size", "15")

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
			// chart.append("text")
		 //        .attr("transform", "rotate(-90)")
		 //        .attr("y", panelMargin.left) // this is actually X direction, because we rotated.
		 //        .attr("x", (-panelYOffset - (panelMargin.top + (panelDims.y / 2))))
		 //        .attr("dy", "-2.7em")
		 //        .style("text-anchor", "middle")
		 //        .text(yLabelText);

	        // Add length label
	        var panelDimsX = panelDims.x * (nucsThisRow / this.nucsPerRow);
			chart.append("text")
		        .attr("x", panelMargin.left + panelDimsX + lengthOffset)
		        .attr("y", panelYOffset + panelMargin.top + panelDims.y)
		        .style("text-anchor", "left")
		        .attr("dy", "1.3em")
		        .text(end);

			// Draw bars
			var barWidth = parseInt(panelDims.x / this.nucsPerRow);
			var bar = chart
				.selectAll("g.nucleotide-measurement-bar r"+rowN)
				.data(dataSlice).enter()
				.append("g")
				.attr("class", "nucleotide-measurement-bar r"+rowN)
				.attr("transform", function(d) { 
					return "translate("+
						(panelMargin.left + xScale(d.position) - (barWidth / 2))+","+
						(panelYOffset + panelMargin.top + yScale(d.measurement))+
					")";
				});

			// Draw the rects for the bars
			bar.append("rect")
				.attr("height", function(d) { 
					return yScale(maxY - d.measurement);
				})
				.attr("width", barWidth)
				.style("fill", "#c33");

			// APPLY STYLES
			var lineStyle = {	
				"fill": "none",
				"stroke": "#000",
				"stroke-width": "1px"
			}
			chart.selectAll("path").style(lineStyle);
	  		chart.selectAll("line").style(lineStyle)




		} // End looping through chart rows
	}
})

/**
 * SearchController handles interactivity for the search module shown on the landing page.
 */
var SearchController = Class.extend({

	// Constructor
	init: function(browserController) {
		this.browserController = browserController;

		$("#search-button").click($.proxy(function(ev) {
 			ev.preventDefault()

 			this.browserController.jumpTo("/search");
			// this.browserController.changeUrl("/search", "Search");
			// $("#transcript-data").empty();
			// $("#search").show()
			// $("#d3nome").hide();
		}, this));

		this.tabController = new TabController([
			new TabControllerTab("search-tab-transcript-id"),
			new TabControllerTab("search-tab-coverage", $.proxy(function() {
				if (this.coverageSearchController == null) {
					this.coverageSearchController = new CoverageSearchController(this.browserController);
				}
			}, this))
		]);

		this.transcriptIDSearchController = new TranscriptIDSearchController(this.browserController);
		this.coverageSearchController = null; // initialises when tab is selected
	},

	show: function() {
		$("#transcript-data").empty();
		$("#search").show()
		$("#d3nome").hide();
		$("#help").hide();
	}
});

var TabController = Class.extend({

	init: function(tabs) {
		this.tabElements = [];
		this.initTabs(tabs);
	},
	initTabs: function(tabs) {
		for (var i = 0; i < tabs.length; i++) {
			if (i == 0) { // first tab is always selected on init
				this.selectedTabID = tabs[i].elementID;
			}
			this.initTab(tabs[i]);
		}
	},
	initTab: function(tab) {
		var element = $("#"+tab.elementID);
		this.tabElements.push(element);
		element.click($.proxy(function(element) {
			var clickedElement = $(element);

			// make sure clicked tab is not same as current tab
			var newSelectedTabID = element.attr("id");
			if (newSelectedTabID == this.selectedTabID) {
				return;
			}
			// keep track of selected tab
			this.selectedTabID = newSelectedTabID;

			// update the tab CSS classes
			for (var i = 0; i < this.tabElements.length; i++) {
				var currElement = this.tabElements[i];
				var currPanelElement = $("#"+currElement.data("ui-id"));

				if (clickedElement.attr("id") != currElement.attr("id")) {
					currElement.removeClass("active")
					currPanelElement.hide()

				} else {
					currElement.addClass("active")
					currPanelElement.show()
				}
			}

			// run tab click callback
			if (tab.clickCallback != null) {
				tab.clickCallback()
			}
		}, this, element))
	}
});

var TabControllerTab = Class.extend({
	init: function(elementID, clickCallback) {
		// clickCallback is an optional field.
		clickCallback = typeof clickCallback !== 'undefined' ? clickCallback : null;

		this.elementID = elementID;
		this.clickCallback = clickCallback;
	}
})


var TranscriptIDSearchController = Class.extend({
	init: function(browserController) {
		this.browserController = browserController;

		var handle = $.proxy(function(ev) {
			if (typeof(ev) !== "undefined") {
				ev.preventDefault();
			}
			var term = $("#search-transcript-id-text").val();
			this.searchTranscriptID(term)
		}, this)

		$("#search-transcript-id-submit").click(handle);
		$('#search-transcript-id-text').on("keypress", $.proxy(function(e) {
        	if (e.keyCode == 13) handle()
        }, this));
	},

	searchTranscriptID: function(term) {
		this.browserController.showLoading();
		$("#search-transcript-id-message").hide();
		$.ajax({
			url: "/ajax/search-transcript/"+term,
			context: this
		}).done(function(results) {
			results = $.parseJSON(results);		
			var exactMatch = false;
			for (var i = 0; i < results.length; i++) {
				if (results[i] == term) {
					// found exact match
					this.browserController.jumpTo("/transcript/"+term);
					// this.browserController.selectTranscript(term);
					return;
				}
			}
			this.browserController.jumpTo("/transcript/"+results[0]);
			// this.browserController.selectTranscript(results[0]);
			this.browserController.hideLoading();

		}).error(function() {
			$("#search-transcript-id-message").html("No transcripts found matching \""+term+"\"")
			$("#search-transcript-id-message").show();
			this.browserController.hideLoading();
		});
	}
});

var CoverageSearchController = Class.extend({
	init: function(browserController) {
		this.browserController = browserController;

		// initialise pagination
		this.browserController.showLoading();
		$.ajax({
			url: "/ajax/get-coverage-page-count", context: this

		}).done(function(pageNum) {
			// insert pagination HTML
			var buf = 
				"	<div id=\"search-coverage-paginator\" class=\"pagination\">" +
                		"<a href=\"#\" class=\"button\" data-action=\"first\">&laquo;</a>" +
                		"<a href=\"#\" class=\"button\" data-action=\"previous\">&lsaquo;</a>" +
                		"<input type=\"text\" class=\"pagination-status\" readonly=\"readonly\" data-max-page=\""+pageNum+"\" />" +
                		"<a href=\"#\" class=\"button\" data-action=\"next\">&rsaquo;</a>" +
                		"<a href=\"#\" class=\"button\" data-action=\"last\">&raquo;</a>" +
                "	</div>" +
				"<div id=\"search-coverage-data\"><!-- filled by paginator AJAX --></div>";
            $("#search-coverage").html(buf);

            // initialise the paginator JS
            $('#search-coverage-paginator').jqPagination({

            	// page change callback
			    paged: $.proxy(function(pageNum) {
			    	this.search(pageNum);
		    	}, this) 
			});

            // retrieve the first page of results.
			this.search(1);
			this.browserController.hideLoading();
		});
	},

	// Grabs transcript coverage data via AJAX and displays it in a div
	search: function(pageNum) {
		this.browserController.showLoading();
		$.ajax({
			url: "/ajax/search-coverage/"+pageNum,
			context: this
		}).done($.proxy(function(results) {
			$("#search-coverage-data").empty();
			$("#search-coverage-data").html(results);

			$(".transcript-row").each($.proxy(function(key, element) {
				var element = $(element)
				element.click($.proxy(function(ev) {
					ev.preventDefault()
					var transcript_id = element.attr("data-transcript-id")
					this.browserController.jumpTo("/transcript/"+transcript_id);
				}, this));
			}, this));
			this.browserController.hideLoading();
		}, this));
	}
});

// Class that handles PCA and structure plotting
var StructureExplorer = Class.extend({
	init: function(browserController) {
		this.browserController = browserController;

		var drawStructureF = $.proxy(function() { this.drawStructure(); }, this)

		this.tabController = new TabController([
			new TabControllerTab("structure-tab-circle-plot", drawStructureF),
			new TabControllerTab("structure-tab-diagram", drawStructureF)
		]);

		this.experimentIDs = [3, 4];
		this.structureData = this.browserController.getJsonFromElement("structure-json")

		this.drawStructurePcas();
		this.initialiseRnaDiagram();
		this.selectedStructure = null;
		this.drawStructure();
	},

	// Set up the forna container - this plots the RNA
	// Also set up button event handlers
	initialiseRnaDiagram: function() {
		if (this.fornaContainer == null) {
			this.fornaContainer = new FornaContainer(
				"#forna-container", {
					'applyForce': true,
					'initialSize': [650, 650]
				}
			);
			this.fornaContainer.setFriction(0.5);
			this.fornaContainer.setCharge(-0.3);
			this.fornaContainer.setGravity(0);
			this.fornaContainer.setPseudoknotStrength(0);
			this.fornaContainer.stopAnimation();
			this.addDmsColours();
		}

		$("#forna-interact-enable").click($.proxy(function(ev) {
			ev.preventDefault();
			this.fornaContainer.startAnimation()
			$("#forna-interact-disable").show()
			$("#forna-interact-enable").hide()
		}, this));

		$("#forna-interact-disable").click($.proxy(function(ev) {
			ev.preventDefault();
			this.fornaContainer.stopAnimation()
			$("#forna-interact-enable").show()
			$("#forna-interact-disable").hide()
		}, this));

		$("#show-dms").click($.proxy(function(ev) {
			ev.preventDefault();
			this.addDmsColours();
		}, this))

		$("#show-ribosome-profiling").click($.proxy(function(ev) {
			ev.preventDefault();
			this.addRibosomeProfilingColours();
		}, this))
	},

	// DMS colours work great
	addDmsColours: function() {
		var measurements = this.getNucleotideMeasurementsFlat(1);

		// manipulate forna into displaying the colours
		this.fornaContainer.addCustomColors({
			color_values: { "": measurements },
			domain: [0, 1],
			range: ["#4f4", "#f44"]
		});
		this.fornaContainer.changeColorScheme("custom");
	},

	// with ribosome, problem is that some colours are really hard to make out.
	// maybe use a log scale to solve this (flatten things out a bit)
	// or a user adjustable threshold scale
	addRibosomeProfilingColours: function() {
		var measurements = this.getNucleotideMeasurementsFlat(2);

		// find max measurement value - use that for max domain
		var max = 0;
		for (var i = 0; i < measurements.length; i++) {
			if (measurements[i] > max) {
				max = measurements[i];
			}
		}

		this.fornaContainer.addCustomColors({
			color_values: { "": measurements },
			domain: [0, max],
			range: ["#fff", "#f00"]
		});
		this.fornaContainer.changeColorScheme("custom");
	},

	getNucleotideMeasurementsFlat: function(experimentID) {
		// get the reactivities as a flat array
		var data = this.browserController
			.getJsonFromElement("nucleotide-measurement-json")[experimentID].data;

		var measurements = [null]; // first element must be ignored
		for (var i = 0; i < data.length; i++) {
			measurements[i + 1] = data[i].measurement;
		}
		return measurements;
	},

	// Get the MFE structure for given experimentID
	getMfe: function(experimentID) {
		var lowestEntry = null;
		var structureData = this.structureData[experimentID].data

		// Find the in vivo structure with the MFE	
		for (var j = 0; j < structureData.length; j++) {
			var currentEntry = structureData[j];
			if (	lowestEntry == null || 
					currentEntry["energy"] < lowestEntry["energy"]) {

				lowestEntry = currentEntry;
			}
		}
		return lowestEntry;
	},

	drawStructurePcas: function() {
		this.drawStructurePca(this.structureData[1], "pca-container-in-silico");
		this.drawStructurePca(this.structureData[2], "pca-container-in-vivo");

		var dlID = window.faGlobals.selectedTid

		new SvgDownloader(
			"pca-container-in-silico-svg", "pca-in-silico-dl", "pca-in-silico_"+dlID+".svg");
		new SvgDownloader(
			"pca-container-in-vivo-svg", "pca-in-vivo-dl", "pca-in-vivo_"+dlID+".svg");

		// attach event handlers to the MFE download buttons
		this.initMfeButton(1, "pca-in-silico-mfe")
		this.initMfeButton(2, "pca-in-vivo-mfe")
	},

	initMfeButton: function(experimentID, buttonID) {
		$("#"+buttonID).click($.proxy(function() {
			this.selectedStructure = this.getMfe(experimentID);
			this.drawStructure();
		}, this));
	},

	// Draws a PCA structure scatter plot
	// http://bl.ocks.org/weiglemc/6185069
	drawStructurePca: function(dataIn, elementID) {
		var svgID = elementID+"-svg";
		var experimentID = dataIn["id"]
		var padding = 0.3; // % margin around the PCA points
		var nTicks = 4;

		// must add container here with button
		var buf = "<svg id=\""+svgID+"\" class=\"structure-pca-chart\"></svg>";

		$("#"+elementID).html(buf)

		dataValues = dataIn["data"];

		var margin = {top: 5, right: 5, bottom: 30, left: 30};
		var totDims = {x: 200, y: 200};
		var panelDims = {
			x: totDims.x - margin.left - margin.right,
			y: totDims.y - margin.left - margin.right
		};

		var energyValue = function(d) { return d.energy; };

		// setup x 
		var xValue = function(d) { return d.pc1; };
		var minX = d3.min(dataValues, xValue);
		var maxX = d3.max(dataValues, xValue);
		var rangeX = maxX - minX;
		var padX = rangeX * padding;
		var xScale = d3.scale.linear()
			.range([0, panelDims.x])
			.domain([
				minX - padX,
				maxX + padX
			]);


		var xMap = function(d) { return xScale(xValue(d)); };
		var xAxis = d3.svg.axis().scale(xScale).ticks(nTicks).orient("bottom");

		// setup y
		var yValue = function(d) { return d.pc2; };
		var minY = d3.min(dataValues, yValue);
		var maxY = d3.max(dataValues, yValue);
		var rangeY = maxY - minY;
		var padY = rangeY * padding;
	    var yScale = d3.scale.linear()
	    	.range([panelDims.y, 0])
	    	.domain([
				minY - padY,
				maxY + padY
			])

	    var yMap = function(d) { return yScale(yValue(d));};
	    var yAxis = d3.svg.axis().scale(yScale).ticks(nTicks).orient("left");

	    // Set up a colour scale
		var numColors = 9;
		var heatmapColour = d3.scale.quantize()
		  	.domain([
		  		d3.min(dataValues, energyValue),
		  		d3.max(dataValues, energyValue),
		  	])
		  	.range(colorbrewer.RdYlGn[numColors]);

		// grab the tooltip element
		var tooltip = d3.select("#structure-pca-chart-tooltip");

		// add the graph canvas to the body of the webpage
		var svg = d3.select("#"+svgID)
		    .attr("width", totDims.x)
		    .attr("height", totDims.y)
		  .append("g")
		    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

		// x-axis
		svg.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate(0," + panelDims.y + ")")
			.call(xAxis)
		  .append("text")
			.attr("class", "label")
			.attr("x", panelDims.x)
			.attr("y", -6)
			.style("text-anchor", "end")
			.text("PC 1");

		// y-axis
		svg.append("g")
			.attr("class", "y axis")
			.call(yAxis)
	      .append("text")
			.attr("class", "label")
			.attr("transform", "rotate(-90)")
			.attr("y", 6)
			.attr("dy", ".71em")
			.style("text-anchor", "end")
			.text("PC 2");

		var lineStyle = {
		  "fill": "none",
		  "stroke": "#000",
		  "shape-rendering": "crispEdges"
		};
		svg.selectAll("path").style(lineStyle)
		svg.selectAll("line").style(lineStyle)

		var showTooltip = function(d) {
			tooltip.transition()
				.duration(0)
				.style("opacity", 1);
			tooltip.html("<i class=\"fa fa-fire\"></i> "+energyValue(d)+" kcal/mol")
				.style("left", (d3.event.pageX) + "px")
				.style("top", (d3.event.pageY) + "px");
		}

		// draw dots
		svg.selectAll(".dot")
			.data(dataValues)
		  .enter().append("circle")
			.attr("class", "dot")
			.attr("r", 5)
			.attr("cx", xMap)
			.attr("cy", yMap)
			.style("fill", function(d) { return heatmapColour(d.energy); })
			.style("stroke", "#000") 
			
			.on("mousemove", showTooltip)
			.on("mouseover", showTooltip)
			.on("mouseout", function(d) {
				tooltip.transition()
					.duration(200)
					.style("opacity", 0);
			})
			.on("click", $.proxy(function(d) {
				this.selectedStructure = d;
				this.drawStructure();
			}, this));
	},

	drawStructure: function() {

		var in_silico_mfe = this.getMfe(1)
		var in_vivo_mfe = this.getMfe(2)

		if (this.selectedStructure == null) {
			this.selectedStructure = in_vivo_mfe;
		}
		window.faGlobals.selectedSid = this.selectedStructure["id"]
		var mfe_txt
		if (this.selectedStructure["id"] == in_vivo_mfe["id"]) {
			mfe_txt = ", <i>in vivo</i> MFE"
		} else if (this.selectedStructure["id"] == in_silico_mfe["id"]) {
			mfe_txt = ", <i>in silico</i> MFE"
		} else {
			mfe_txt = ""
		}


		$("#forna-energy").html(this.selectedStructure["energy"]+" kcal/mol"+mfe_txt);

		if (this.tabController.selectedTabID == "structure-tab-diagram") {
			this.drawStructureDiagram();
			$("#circle-plot-legend").hide()
		} else {
			this.drawCirclePlot();
			$("#circle-plot-legend").show()
		}
	},

	drawStructureDiagram: function() {
		var structureID = this.selectedStructure["id"];

		// this should be done in the constructor really.
		// or have a getStructure method that attached to the object with a callback
		this.browserController.showLoading();
		$.ajax({
			url: "/ajax/structure-diagram/"+structureID, 
			context: this
		}).done(function(data) {

			// This data includes sequence string, dot bracket structure
			// and 2d diagram coords.
			data = JSON.parse(data);

			g = new RNAGraph(data["sequence"], data["structure"], "rna")
		        .elementsToJson()
		        .addPositions('nucleotide', data["coords"])
		        .addLabels(1) // 1 = start
		        .reinforceStems()
		        .reinforceLoops()
		        .connectFakeNodes();

		    this.fornaContainer.clearNodes();// remove the previous diagram
			this.fornaContainer.addRNAJSON(g, true); // generate new diagram
			this.browserController.hideLoading();
			var dlID = window.faGlobals.selectedTid+"-"+window.faGlobals.selectedSid;
			new SvgDownloader("plotting-area", "fornac-dl-button", "structure-diagram_"+dlID+".svg");
		});
	},

	drawCirclePlot: function() {
		var structureID = this.selectedStructure["id"];
		// this should be done in the constructor really.
		// or have a getStructure method that attached to the object with a callback
		this.browserController.showLoading();
		$.ajax({
			url: "/ajax/structure-circle-plot/"+structureID, 
			context: this
		}).done(function(data) {

			$("#circle-plot").empty();
			$("#circle-plot").append(
				"<a href=\"javascript:void(0)\" id=\"circle-plot-dl-button\" class=\"button svg\">"+
					"<i class=\"fa fa-file-image-o\"></i>"+
				"</a>"
			)
			data = JSON.parse(data);

			// find min and max bpp values
			var min_bpp = null;
			for (var i = 0; i < data.length; i++) {
				var bpp = data[i].bpp
				if (min_bpp == null) min_bpp = bpp;
				if (bpp == null) continue;
				if (bpp < min_bpp) min_bpp = bpp;
			}
			var max_bpp = 0

			// now add the circle plot legend
			$("#bpp_low").html(Math.round(min_bpp * 10) / 10)
			$("#bpp_high").html("0.0")

			var nTicks = 10;
			var nNucs = data.length;
			var getColour = function(position) {
				var bpp = data[position].bpp;
				if (bpp == null) {
					return d3.hsl(hue, 0, 0.5);
				}
				var diff = max_bpp - min_bpp;
				var intensity = (bpp - min_bpp) / diff

				// see docs: https://github.com/mbostock/d3/wiki/Colors
				var hue =  intensity * 135; // 360;
				return d3.hsl(hue, 0.8, 0.35);
			}

			var svgDims = 630,
				diameter = 550,
			    radius = diameter / 2,
			    innerRadius = radius - 18;

			var offset = (svgDims - diameter) / 2

			var cluster = d3.layout.cluster()
			    .size([360, innerRadius])
			    .sort(null)
			    .value(function(d) { return d.size; });

			var bundle = d3.layout.bundle();

			var line = d3.svg.line.radial()
			    .interpolate("bundle")
			    .tension(0.15)
			    .radius(function(d) { return d.y; })
			    .angle(function(d) { return d.x / 180 * Math.PI; });

			var svg = d3.select("#circle-plot").append("svg")
				.attr("id", "circle-plot-svg")
			    .attr("width", svgDims)
			    .attr("height", svgDims)
			  .append("g")
			    .attr("transform", "translate(" + (radius + offset) + "," + (radius + offset) + ")");

			// svg.append("circle")
			// 	.attr("class", "circleplot-circle")
   //  			.attr("r", radius);

			svg.append("circle")
				.attr("class", "circleplot-circle")
				.attr("stroke", "#000")
				.attr("fill", "#fff")
    			.attr("r", innerRadius);

			var link = svg.append("g").selectAll(".circleplot-link"),
			    node = svg.append("g").selectAll(".circleplot-node");


			var nodes = cluster.nodes(prepareData(data)),
				links = prepareLinks(nodes);

			link = link
				.data(bundle(links))
				.enter().append("path")
				.each(function(d) { d.source = d[0], d.target = d[d.length - 1]; })
				.attr("class", "circleplot-link")
				.attr("stroke", "#000")
				// .attr("stroke-opacity", ".8")

				.attr("fill", "none")
				.attr("d", line)
				.style("stroke", function(d) {
					return getColour(d[0].key);
				});

			var tickWidth = Math.floor(data.length / nTicks);

			// Calculate ticks by building a throwaway axis
			var axis = d3.svg.axis()
				.scale(d3.scale.linear().domain([0, nNucs]))
				.ticks(10);

			// Retrieve tick values to compare
			var tickValues = axis.scale().ticks(axis.ticks()[0])

			// Build some tick node groups and move them to the right spots, with the 
			// correct rotations and translations
			ticks = node
				.data(nodes.filter(function(n) { 
					// this ensures ticks start at 1 and then goes up in (e.g.) 
					// 100 nuc intervals
					var key = n.key + 1
					if (key == 1) key = 0;
					return !n.children && tickValues.indexOf(key) > -1;
				}))
				.enter().append("g")
				.attr("transform", function(d) { 
					return 	"rotate(" + (d.x - 90) + ")"+
							"translate(" + (d.y) + ",0)"+
							"rotate(90)"
				})

			// Attach a text label to each node
			ticks.append("text")
				.attr("class", "circleplot-node")
				.attr("dy", ".31em")
				.attr("transform", function(d) { return "translate(0, -15)"; })
				.style("text-anchor", function(d) {
					return "middle"
				})

				.text(function(d) { return d.key + 1; })

			// Attach a tick line to each tick node
			ticks.append("path")
				.attr("class", "circleplot-tick")
				.attr("stroke", "#000")
				.attr("d", "M 0 -6 L 0 6")

			this.browserController.hideLoading();

			var dlID = window.faGlobals.selectedTid+"-"+window.faGlobals.selectedSid;
			new SvgDownloader("circle-plot-svg", "circle-plot-dl-button", "circle-plot_"+dlID+".svg");

			// reorganise the data a bit
			function prepareData(rawData) {
				map = {
					name: "",
					children: [],
				};
				for (var i = 0; i < rawData.length; i++) {
					map.children[i] = {
						key: i, // TODO is this needed?
						link: rawData[i].link,
						name: i,
						label: rawData[i].label,
						parent: map
					}
				}

				return map;
			}

			// make links between nodes
			function prepareLinks(nodes) {
				var map = {},
				links = [];
				nodes.forEach(function(d) {
					map[d.name] = d;
				});
				nodes.forEach(function(d) {
					if (d.link != null) {
						links.push({
							source: map[d.name], 
							target: map[d.link]
						});
					};
				});
				return links;
			}
		});
	}
});

// Binds a button to SVG download
// see https://stackoverflow.com/questions/23218174/how-do-i-save-export-an-svg-file-after-creating-an-svg-with-d3-js-ie-safari-an
var SvgDownloader = Class.extend({

	// Connect the download link to an SVG download option
	init: function(svgID, linkID, filename) { // add filename
		$("#"+linkID).click($.proxy(function() {

			//get svg element.
			var svg = document.getElementById(svgID);

			//get svg source.
			var serializer = new XMLSerializer();
			var source = serializer.serializeToString(svg);

			//add name spaces.
			if(!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)){
			    source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
			}
			if(!source.match(/^<svg[^>]+"http\:\/\/www\.w3\.org\/1999\/xlink"/)){
			    source = source.replace(/^<svg/, '<svg xmlns:xlink="http://www.w3.org/1999/xlink"');
			}

			//add xml declaration
			source = '<?xml version="1.0" standalone="no"?>\r\n' + source;

			// download the file
			var blob = new Blob([source], {type: "image/svg+xml;charset=utf-8"});
		    saveAs(blob, filename);
		}, this));
	}
});




