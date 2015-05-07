// Ensembl REST API Gene model
Genoverse.Track.Model.Gene.DBGene = Genoverse.Track.Model.Gene.extend({

  url: "http://rnabrowser.dev/ajax/genome-browser/genes?chr=__CHR__&start=__START__&end=__END__",
  
  // The url above responds in json format, data is an array
  // We assume that parents always preceed children in data array, gene -> transcript -> exon
  // See rest.ensembl.org/documentation/info/feature_region for more details
  parseData: function (data) {
    for (var i = 0; i < data.length; i++) {
      var feature = data[i];
      
      if (feature.feature_type === 'gene' && !this.featuresById[feature.id]) {
        feature.label       = feature.external_name || feature.id;
        feature.transcripts = [];
        
        this.insertFeature(feature);
      }
    }
  }
});
