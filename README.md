Map.setOptions('satellite')
Map.centerObject(geometry,13)


// set up filter Landsat composite

function maskClouds(image) {
  
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
    var cloudShadowBitMask = ee.Number(2).pow(3).int();
    var cloudsBitMask = ee.Number(2).pow(5).int();  
    
    // Get the pixel QA band.
    var qa = image.select('pixel_qa');
    
     // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).and(qa.bitwiseAnd(cloudsBitMask).eq(0)); 
  
  // Return the masked image, scaled to [0, 1].
  return image.updateMask(mask).divide(10000).copyProperties(image, ["system:time_start"]);
}

//2.2) Adding Spectral Indices
///////////////////////////////

// This function maps spectral indices for Mangrove Mapping using Landsat 8 Imagery
var addIndicesL8 = function(img) {
  // NDVI
  var ndvi = img.normalizedDifference(['B5','B4']).rename('NDVI');
  // NDMI (Normalized Difference Mangrove Index - Shi et al 2016 - New spectral metrics for mangrove forest identification)
  var ndmi = img.normalizedDifference(['B7','B3']).rename('NDMI');
  // MNDWI (Modified Normalized Difference Water Index - Hanqiu Xu, 2006)
  var mndwi = img.normalizedDifference(['B3','B6']).rename('MNDWI');
  // SR (Simple Ratio)
  var sr = img.select('B5').divide(img.select('B4')).rename('SR');
  // Band Ratio 54
  var ratio54 = img.select('B6').divide(img.select('B5')).rename('R54');
  // Band Ratio 35
  var ratio35 = img.select('B4').divide(img.select('B6')).rename('R35');
  // GCVI
  var gcvi = img.expression('(NIR/GREEN)-1',{
    'NIR':img.select('B5'),
    'GREEN':img.select('B3')
  }).rename('GCVI');
  return img
    .addBands(ndvi)
    .addBands(ndmi)
    .addBands(mndwi)
    .addBands(sr)
    .addBands(ratio54)
    .addBands(ratio35)
    .addBands(gcvi);
};

//2.3) Filter Landsat data by Date and Region
/////////////////////////////////////////////

// Temporal Parameters

// Select the desired central year here
var year = 2013; 

// Start date will be set one year before the central year
var startDate = (year-1)+'-01-01'; 

// End date will be set to one year later than the central year.
var endDate = (year+1)+'-12-31'; 

//2.4) Apply filters and masks to Landsat 8 imagery
////////////////////////////////////////////////////

var l8 = L8.filterDate(startDate,endDate)
// Mask for clouds and cloud shadows
    .map(maskClouds)
//Add the indices
    .map(addIndicesL8)
    
//2.5) Composite the Landsat image collection
/////////////////////////////////////////////

//You can composite on a per pixel, per-band basis using .median()
// OR with quality bands like .qualityMosaic('NDVI')

var composite = l8
              // Uses the median reducer
              .median() 
              // Clips the composite to our area of interest
              .clip(geometry); 

//2.6) Mask to areas of low elevation and high NDVI and MNDWI
/////////////////////////////////////////////////////////////

// Clip SRTM data to region
var srtmClip = SRTM.clip(geometry);

//Mask to elevations less than 65 meters
var elevationMask = srtmClip.lt(65);

//Used the NDVI and MNDWI bands to create masks
var NDVIMask = composite.select('NDVI').gt(0.25);
var MNDWIMask = composite.select('MNDWI').gt(-0.50);

//Apply the masks
var compositeNew = composite
                        .updateMask(NDVIMask)
                        .updateMask(MNDWIMask)
                        .updateMask(elevationMask)
                        
//2.7) Display results
///////////////////////

//Select bands and parameters for visualization
var visPar = {bands:['B5','B6','B4'], min: 0, max: 0.35}; 

//Add layer to map
Map.addLayer(compositeNew.clip(geometry), visPar, 'Landsat Composite 2013')

///////////////////////////////////////////////////////////////
//          3) Construct Random Forest Model                 //
///////////////////////////////////////////////////////////////

//3.1) Prepare training data and predictors
////////////////////////////////////////////

//After drawing training polygons, merge them together
var classes = Mangrove.merge(NonMangrove)

//Define the bands you want to include in the model
var bands = ['B5','B6','B4','NDVI','MNDWI','SR','GCVI']

//Create a variable called image to select the bands of interest and clip to geometry
var image = compositeNew.select(bands).clip(geometry)
   
//Assemble samples for the model
var samples = image.sampleRegions({
    collection: classes, // Set of geometries selected for training
    properties: ['landcover'], // Label from each geometry
    scale: 30 // Make each sample the same size as Landsat pixel
    }).randomColumn('random'); // creates a column with random numbers
    
//Here we randomly split our samples to set some aside for testing our model's accuracy
// using the "random" column we created
var split = 0.8; // Roughly 80% for training, 20% for testing.
var training = samples.filter(ee.Filter.lt('random', split)); //Subset training data
var testing = samples.filter(ee.Filter.gte('random', split)); //Subset testing data


//Print these variables to see how much training and testing data you are using
    print('Samples n =', samples.aggregate_count('.all'));
    print('Training n =', training.aggregate_count('.all'));
    print('Testing n =', testing.aggregate_count('.all'));

//3.2) Begin Random Forest Classification
/////////////////////////////////////////

//.smileRandomForest is used to run the model. Here we run the model using 100 trees
// and 5 randomly selected predictors per split ("(100,5)")
    var classifier = ee.Classifier.smileRandomForest(100,5).train({ 
    features: training.select(['B5','B6','B4','NDVI','MNDWI','SR','GCVI', 'landcover']), //Train using bands and landcover property
    classProperty: 'landcover', //Pull the landcover property from classes
    inputProperties: bands
    });

//3.3) Test the accuracy of the model
//////////////////////////////////////

    var validation = testing.classify(classifier);
    var testAccuracy = validation.errorMatrix('landcover', 'classification');
    print('Validation error matrix RF: ', testAccuracy);
    print('Validation overall accuracy RF: ', testAccuracy.accuracy());

//3.4) Classify the Landsat composite using the Random Forest model
///////////////////////////////////////////////////////////////////

    var classifiedrf = image.select(bands) // select the predictors
                      .classify(classifier); // .classify applies the Random Forest
                      
//The model results may be "noisy". To reduce noise, create a mask to mask
// unconnected pixels
    var pixelcount = classifiedrf.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
    var countmask = pixelcount.select(0).gt(25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
    var classMask = classifiedrf.select('classification').gt(0)
    var classed= classifiedrf.updateMask(countmask).updateMask(classMask)

//3.5) Map results
////////////////////

//Add classification to map
Map.addLayer (classed, {min: 1, max: 1, palette:'blue'}, 'Mangrove 2013');

//Use reduceRegion with a Sum reducer to calculate total area
var get2013 = classed.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:geometry,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
print(get2013, 'Mangrove Extent 2013 in ha')


//Overall Accuracy Mangrove Forest
var uji_akurasi = uji_akurasi_Mangrove2013.merge(uji_akurasi_NonMangrove2013)

var validasi = classed.sampleRegions({
  collection: uji_akurasi,
  properties: ['Landcover'],
  scale: 30,
});
print(validasi);

var akurasi = validasi.errorMatrix('Landcover', 'classification');
print('Confussion Matrix', akurasi);
print('Overall Accuracy', akurasi.accuracy());

Export.image.toDrive({
  image: classed,
  description: '2013GiliMenoMangroveExtent',
  region: geometry,
  scale: 30,
  maxPixels: 1e13
  });
