var express = require('express');
var odex = require("odex");

var router = express.Router();


/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('index', { title: 'Methylmercury Burden Estimator' });
});

//router.get('/draw', function(req, res, next) {
//  res.json({ message: 'Hello World' });
//});



var OrganScaleFactor = function(weight, meanWeight) {
  return weight/meanWeight;
}

var FlowScaleFactor = function(weight, meanWeight) {
  return (weight**0.75)/(meanWeight**0.75);
}

var RateScaleFactor = function(weight, meanWeight) {
  return (weight**(-0.25))/(meanWeight**(-0.25));
}

var slowDoses = function(t, period, dose, halflife, end) {
  if(t>end)
    return 0;
  var time = t % period;
  return dose * (2.7182818284**(-time/halflife) / halflife);
}

var param_setup = function(patientWeight, currentTarget) {
  var weightArray = [19.28, 83.12, 68.91];
  var meanWeight = weightArray[currentTarget];

  var totalWellPerfusedVolArray = [3.26, 8.18, 6.85];
  var totalPoorlyPerfusedVolArray = [5.6, 36, 25];
  var brainVolArray = [1.39, 1.34, 1.20];
  var kidneyVolArray = [0.09, 0.32, 0.26];
  var fatVolArray = [10.40, 32.19, 34.78];
   
  var liverVolArray = [0.50, 1.57, 1.35];
  var gutVolArray = [0.28, 1.23, 1.17];
  var lumenVolArray = [0.42, 0.9, 0.83]; 
  var plasmaVolArray =  [2.4, 5.3, 3.9];
  var bloodcellVolArray = [2.4, 5.3, 3.9];

  var skeletalMuscleVolArray = [4.9, 32, 21];

  var KFeArray = [0.00292, 0.00625, 0.005 ]
  var totalWellPerfusedFlowArray = [2.6, 8, 7.44];
  var totalPoorlyPerfusedFlowArray = [2.69, 3.79, 3.77];
  var brainFlowArray = [0.72, 0.68, 0.63];
  var kidneyFlowArray = [0.33, 1.17, 0.85];
  var liverFlowArray = [0.46, 1.32, 1.35];
  var gutFlowArray = [0.22, 0.93, 0.91];
  var skeletalMuscleFlowArray = [0.15, 0.95, 0.63];
  var fatFlowArray = [0.28, 0.68, 1.04];

  var totalWellPerfusedVol = totalWellPerfusedVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);
  var totalPoorlyPerfusedVol = totalPoorlyPerfusedVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);
  var VBr = brainVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);
  var VKi = kidneyVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);
  var VLv = liverVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);
  var VGt = gutVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);
  var VGL = lumenVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);
  var VPl = plasmaVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight) * 0.6;
  var VMu = skeletalMuscleVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);
  var VRbc = bloodcellVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight) * 0.4;
  var VFa = fatVolArray[currentTarget] * OrganScaleFactor(patientWeight, meanWeight);

  var totalWellPerfusedFlow = totalWellPerfusedFlowArray[currentTarget] * FlowScaleFactor(patientWeight, meanWeight) * 0.6 * 60;
  var totalPoorlyPerfusedFlow = totalPoorlyPerfusedFlowArray[currentTarget] * FlowScaleFactor(patientWeight, meanWeight) * 0.6 * 60;
  var QBr = brainFlowArray[currentTarget] * FlowScaleFactor(patientWeight, meanWeight) * 0.6 * 60;
  var QKi = kidneyFlowArray[currentTarget] * FlowScaleFactor(patientWeight, meanWeight) * 0.6 * 60;
  var QFa = fatFlowArray[currentTarget] * FlowScaleFactor(patientWeight, meanWeight) * 0.6 * 60;
  var QLv = liverFlowArray[currentTarget] * FlowScaleFactor(patientWeight, meanWeight) * 0.6 * 60;
  var QGt = gutFlowArray[currentTarget] * FlowScaleFactor(patientWeight, meanWeight) * 0.6 * 60;
  var QMu = skeletalMuscleFlowArray[currentTarget] * FlowScaleFactor(patientWeight, meanWeight) * 0.6 * 60;
  var VSp = totalPoorlyPerfusedVol - VMu;
  var VRp = totalWellPerfusedVol - VBr - VLv - VKi - VGt;
  var QSp = (totalPoorlyPerfusedFlow - QMu);
  var QRp = (totalWellPerfusedFlow - QBr - QLv - QKi  - QGt);
  var kRPl = 0.3 * RateScaleFactor(patientWeight, meanWeight);
  var kPlR = 3 * RateScaleFactor(patientWeight, meanWeight);
            //  0      1    2    3    4    5    6    7    8     9           10
  return [meanWeight, VBr, VKi, VLv, VGt, VGL, VPl, VMu, VRbc, VFa, totalWellPerfusedFlow, 

       //     11             12   13   14   15   16   17   18   19   20   21   22     23
    totalPoorlyPerfusedFlow, QBr, QKi, QFa, QLv, QGt, QMu, VSp, VRp, QSp, QRp, kRPl, kPlR];
}

var SolveSystem = function(KFe, KBi, kGLI, kEx, kLvI, kAbs, weight, btype, fishPpm, fishVol, dose) {

  //Set up params
  var PBr = 12 * 3;
  var PMu = 12 * 2; 
  var PRp = 12 * 1;
  var PSp = 12 * 2;
  var PLv = 12 * 5;
  var PGt = 12 * 1;
  var PKi = 12 * 4;
  var PFa = 10 * 0.15;

  //console.log("import: " + KFe, KBi, kGLI, kEx, kLvI, kAbs, btype);

  var params = param_setup(Number(weight), btype - 1);

  //Set up system
                      //  0      1    2    3    4    5    6    7    8     9           10
  var derivs = function(weight, VBr, VKi, VLv, VGt, VGL, VPl, VMu, VRbc, VFa, totalWellPerfusedFlow, 

                           //     11             12   13   14   15   16   17   18   19   20   21   22     23
                        totalPoorlyPerfusedFlow, QBr, QKi, QFa, QLv, QGt, QMu, VSp, VRp, QSp, QRp, kRPl, kPlR,

                      // 24  25   26    27    28   29
                        KFe, KBi, kGLI, kEx, kLvI, kAbs, end, dose) {

    //console.log(//VBr,VKi,VLv,VGt,VGL, VPl, VMu, VRbc, VFa, totalWellPerfusedFlow, 

      //totalPoorlyPerfusedFlow, QBr, QKi, QFa, QLv, QGt, QMu, VSp, VRp, QSp, QRp, kRPl, kPlR, 
      //KFe, KBi, kGLI, kEx, kLvI, kAbs);

      return function(x, y) {
      return[
      
        //APl'[t]   0
        (QLv * y[6]/(VLv * PLv) - (QLv - QGt) * y[0]/VPl) + 
          QBr * (y[1]/(VBr * PBr) - y[0]/VPl) + 
          QSp * (y[5]/(VSp * PSp) - y[0]/VPl) + 
          QRp * (y[4]/(VRp * PRp) - y[0]/VPl) + 
          QKi * (y[9]/(VKi * PKi) - y[0]/VPl) + 
          QMu * (y[2]/(VMu * PMu) - y[0]/VPl) +  
          QFa * (y[3]/(VFa *PFa) - y[0]/VPl) - QGt * y[0]/VPl + 
          kRPl * (y[10]) - kPlR * (y[0]) + 
          slowDoses(x, 24*7, dose, 3, end),

        //ABr'[t]   1
        -QBr * (y[1]/(VBr * PBr) - y[0]/VPl),


        //AMu'[t]   2 
        -QMu * (y[2]/(VMu * PMu) - y[0]/VPl),


        //AFa'[t]   3
        -QFa * (y[3]/(VFa * PFa) - y[0]/VPl),


        //ARp'[t]   4
        -QRp * (y[4]/(VRp * PRp) - y[0]/VPl),


        //ASp'[t]   5
        -QSp * (y[5]/(VSp * PSp) - y[0]/VPl),


        //ALv'[t]   6
        -(QLv * y[6]/(VLv * PLv) - (QLv - QGt) * y[0]/VPl - 
            QGt * y[7]/(VGt * PGt)) - KBi * y[0]/VPl - kLvI * y[6],


        //AGt'[t]   7
        -QGt * (y[7]/(VGt * PGt) - y[0]/VPl) + 
          kAbs * y[8] - kEx * y[7],


        //AGL'[t]   8 
        KBi * y[0]/VPl - kAbs * y[8] - KFe * y[8]/VGL - kGLI * y[8] + 
          kEx * y[7],

        //AKi'[t]   9
        -QKi * (y[9]/(VKi * PKi) - y[0]/VPl),

        //ARBC'[t]  10
        -kRPl * y[10] + kPlR * y[0],

        //ABiTransTot'[t] 11
        KBi * y[0]/VPl ,

        //AGtTransTot'[t] 12
        kEx * y[7]
      ];

    };
  };
  var s3 = new odex.Solver(13);
  s3.absoluteTolerance = s3.relativeTolerance = 1e-8;
  s3.denseOutput = true;
  s3.maxSteps = 300000;

  var times = [];
  var blood_cons = [];
  var brain_cons = [];
  var body_cons = [];


  var n_steps = 3*800;
  var end3 = 1400;
  var sample_dist = 3.2;
  /*
  s1.solve(derivs(params[0], params[1], params[2], params[3], params[4], params[5], params[6], 
    params[7], params[8], params[9], params[10], params[11], params[12], params[13], params[14],params[15], 
    params[16], params[17], params[18], params[19], params[20], params[21], params[22], params[23],
    KFe, KBi, kGLI, kEx, kLvI, kAbs, end1, dose), 
    0,
    [0,0,0,0,0,0,0,0,0,0,0,0,0], 
    n_points,
    s1.grid(sample_dist, function(x,y) {
      plasma_cons1.push(y[0] / params[6]);
      rbc_cons1.push(y[10] / params[8]);
  }));

  s2.solve(derivs(params[0], params[1], params[2], params[3], params[4], params[5], params[6], 
    params[7], params[8], params[9], params[10], params[11], params[12], params[13], params[14],params[15], 
    params[16], params[17], params[18], params[19], params[20], params[21], params[22], params[23],
    KFe, KBi, kGLI, kEx, kLvI, kAbs, end2, dose), 
    0,
    [0,0,0,0,0,0,0,0,0,0,0,0,0], 
    n_points,
    s2.grid(sample_dist, function(x,y) {
      plasma_cons2.push(y[0] / params[6]);
      rbc_cons2.push(y[10] / params[8]);
  }));
  */
  s3.solve(derivs(params[0], params[1], params[2], params[3], params[4], params[5], params[6], 
    params[7], params[8], params[9], params[10], params[11], params[12], params[13], params[14],params[15], 
    params[16], params[17], params[18], params[19], params[20], params[21], params[22], params[23],
    KFe, KBi, kGLI, kEx, kLvI, kAbs, end3, dose), 
    0,
    [0,0,0,0,0,0,0,0,0,0,0,0,0], 
    n_steps,
    s3.grid(sample_dist, function(x,y) {
      blood_cons.push((y[0] + y[10]) / (params[6] + params[8])); //units are micrograms per liter.
      brain_cons.push(y[1] / params[1]);
      body_cons.push(y[2] / params[7]);
      times.push(x);
  }));

  var p1 = Math.floor(n_steps/sample_dist - 100);
  var p2 = Math.floor(n_steps/sample_dist) - 1;
  var max_blood_con = Math.max.apply(Math, blood_cons);

  var halflife = (((p2 - p1) * sample_dist) / (Math.log(blood_cons[p1] / blood_cons[p2]) / Math.log(2))) / 24;
  var equ_time = halflife * Math.log((dose/(24*7)) / (0.05 * max_blood_con / (24*7))) / Math.log(2);

  return [times, blood_cons, brain_cons, body_cons, halflife, equ_time];

}











router.post('/draw', function (req, res) {
  var KFe = 0.00661087;
  var KBi = 1.71883;
  var kGLI = 0.0136666;
  var kEx = 0.273333;
  var kLvI = 0.000911109;
  var kAbs = 0.273333;

  var weight = Number(req.body.weight) * Number(req.body.weightUnits);
  var btype = Number(req.body.btype);
  var fishPpm = Number(req.body.ftype);
  var fishVol = Number(req.body.meal_size) * Number(req.body.mealUnits);
  if(isNaN(weight) || isNaN(btype) || isNaN(fishPpm) || isNaN(fishVol)) {
    res.render('non-numeric-input');
  }
  else {
    var dose = fishPpm * fishVol; //units are micrograms.

    var sysSolution = SolveSystem(KFe, KBi, kGLI, kEx, kLvI, kAbs, weight, btype, fishPpm, fishVol, dose);

    res.render('draw', { data: 'Testing...', timeVals: sysSolution[0], blood_vals: sysSolution[1], brain_vals: sysSolution[2], 
    body_vals: sysSolution[3],  halflife: sysSolution[4], equ_time: sysSolution[5] });
  }
})

router.post('/advanced', function (req, res) {

  res.render('advanced');
})

router.post('/drawAdvanced', function (req, res) {
  var KFe = Number(req.body.KFe);
  var KBi = Number(req.body.KBi);
  var kGLI = Number(req.body.kGLI);
  var kEx = Number(req.body.kEx);
  var kLvI = Number(req.body.kLvI);
  var kAbs = Number(req.body.kAbs);

  var weight = Number(req.body.weight) * Number(req.body.weightUnits);
  var btype = Number(req.body.btype);
  var fishPpm = Number(req.body.ftype);
  var fishVol = Number(req.body.meal_size) * Number(req.body.mealUnits);
  var dose = fishPpm * fishVol;

  if(isNaN(weight) || isNaN(btype) || isNaN(fishPpm) || isNaN(fishVol) || isNaN(KFe) || isNaN(KBi)
   || isNaN(kGLI) || isNaN(kEx) || isNaN(kLvI) || isNaN(kAbs)) {
    res.render('non-numeric-input');
  }
  else {
    console.log(KFe, KBi, kGLI, kEx, kLvI, kAbs, weight, btype, fishPpm, fishVol, dose);
    var sysSolution = SolveSystem(KFe, KBi, kGLI, kEx, kLvI, kAbs, weight, btype, fishPpm, fishVol, dose);
    
    res.render('drawAdvanced', { data: 'Testing...', timeVals: sysSolution[0], blood_vals: sysSolution[1], brain_vals: sysSolution[2], 
    body_vals: sysSolution[3],  halflife: sysSolution[4], equ_time: sysSolution[5] });
  }
})

module.exports = router;
