// Custom gene loading component
// Talks to RNA browser endpoint to fetch genome data
// @author Matthew Norris

// TODO rename to GffSQL
Genoverse.Track.DBGene = Genoverse.Track.extend({
  id     : 'genes',
  name   : 'Genes',
  height : 200,

  populateMenu: function (feature) {
    // get the transcript ID

    // do something with the transcript ID
    console.log(feature)
  },

  // 2000000: { // This one applies when > 2M base-pairs per screen
  //   labels : false
  // },
  // 100000: {
  //   labels : false,
  //   model  : Genoverse.Track.Model.Gene.Ensembl,
  //   view   : Genoverse.Track.View.Gene.Ensembl
  // },
  // 1: { // > 1 base-pair, but less then 100K
  //   labels : true,
  //   model  : Genoverse.Track.Model.Transcript.GFF3,
  //   view   : Genoverse.Track.View.Transcript.Ensembl
  // }

  2000000: { // This one applies when > 2M base-pairs per screen
    labels : true,
    model  : Genoverse.Track.Model.Transcript.DBTranscript,
    view   : Genoverse.Track.View.Transcript.Ensembl
  },
  100000: {
    labels : true,
    model  : Genoverse.Track.Model.Transcript.DBTranscript,
    view   : Genoverse.Track.View.Transcript.Ensembl
  },
  1: { // > 1 base-pair, but less then 100K
    labels : true,
    model  : Genoverse.Track.Model.Transcript.DBTranscript,
    view   : Genoverse.Track.View.Transcript.Ensembl
  }
});
