/**
 * Macro for data evaluation of DIRK_P700_el
 * data on PhotosynQ.org
 * by: Andre Potthoff
 * created: 2025-06-07
 * based on macro from David M. Kramer (RIDES2.0)
 */

// Define the output object here
var output = {}; 

// ##########################################
// ###                                    ###
// ###           P700_analysis            ###
// ###                                    ###
// ##########################################

// ASSUMPTIONS about the protocol:
    // "v_arrays": [
    //   [
    //     1000,  <-- the time distance between the measurement pulses in µs (constant - I would have used save_trace_time_scale but it does not work)
    //     150,   <-- the number of pulses in the dark period
    //     30     <-- the number of pulses in the dark period. determines the length of the dark-period.
    //   ],
    //   [
    //     50,    <-- all the PARs tested (sequentially after another)
    //     100,
    //     500,
    //     1000,
    //     5000,
    //     10000,
    //     50000
    //   ]
    // ],


var DIRK_P700_trace_beginning = json.v_arrays[0][1]; // initial acclimation time in seconds (yes we have pre-illumination but still)
var DIRK_P700_length_of_baseline = json.v_arrays[0][1]; // first 100 pulses under light
var DIRK_P700_num_pulses_dark_period = json.v_arrays[0][2];// length of relaxation kinetics
var DIRK_P700_recovery = json.v_arrays[0][1]; // 100 pulses recovery, typically curves go back to baseline much faster
var DIRK_P700_number_of_subtraces = 6; // baseline, dark, and recovery are repeated 7 times
var DIRK_P700_end_of_dark_period = DIRK_P700_length_of_baseline + DIRK_P700_num_pulses_dark_period; // end of the dark period in number of pulses
var DIRK_P700_length_of_subtrace = DIRK_P700_length_of_baseline + DIRK_P700_num_pulses_dark_period + DIRK_P700_recovery; 
var DIRK_P700_total_number_of_measurement_pulses = DIRK_P700_trace_beginning + DIRK_P700_number_of_subtraces * DIRK_P700_length_of_subtrace; // total length of the trace in number of measurement pulses
var DIRK_P700_initial_dark_period_pulses_used_for_fitting = 30;

var DIRK_P700_time_beween_pulses = json.v_arrays[0][0] / 1000; // get the time axis from the json file in ms (from µs)
var DIRK_P700_PAR_levels_tested = json.v_arrays[1];
var DIRK_P700_number_of_tested_PARS = DIRK_P700_PAR_levels_tested.length;
var DIRK_P700_only_one_PAR_level_tested = DIRK_P700_number_of_tested_PARS == 1;

var DIRK_P700_a_initial = 0.2;
var DIRK_P700_b_initial = -0.002;    // the choice of these initial values is very critical for relaible fitting. (writing this in pain). these values have shown to work well
var DIRK_P700_c_initial = 1500;
var DIRK_P700_fitted_curve=[];
var DIRK_P700_out_P700_amplitude = [];
var DIRK_P700_out_k_PSI = [];               //  rate constant of PSI P700 re-reduction
var DIRK_P700_out_P700_tau = [];            // 1 / k_PSI (unit: s)
var DIRK_P700_out_v_PSI = [];               // intial proton flux through ATP Synthase when dark_period starts (unit: mAU s^-1)
var DIRK_P700_out_fit_R2 = [];              // R^2 value of the fit
var DIRK_P700_out_raw_transmission = [];
var DIRK_P700_out_average_transmission = [];
var DIRK_P700_out_averaged_absorbtion = [];
var DIRK_P700_out_relaxation_fit_data_and_fit = [];
var DIRK_P700_out_residuals_of_exponential_fit = [];
var DIRK_P700_out_dark_relaxation_window = [];
var DIRK_P700_out_fitted_curve_exp_decay = [];
var DIRK_P700_out_b = [];


var decay = function(x,a,b,c){
  return b + a * Math.exp( -x/c);
};

function unwrap_if_single(input, stringify = false) {
    if (Array.isArray(input) && input.length === 1) {
        out = input[0];
        if (stringify){
          return out.join(" ");
        }
        return out;
    }
    if (stringify){
      return input.join(" ");
    }
    return input;
}

//////////////////////////////////////////////
// GET the important objects from the input //
//////////////////////////////////////////////

var DIRK_P700_traces = GetProtocolByLabel("DIRK_P700_hl", json, true); // get all traces from json

// create a time axis for the relaxation data, containing 21 values from 0 to "dark_period"
var DIRK_P700_time_axis = [];
for (i = 0; i < DIRK_P700_initial_dark_period_pulses_used_for_fitting+1; i++){
  DIRK_P700_time_axis.push(i*DIRK_P700_time_beween_pulses);
}

////////////////////////////////////////////
// Loop through all iterations of the set //
////////////////////////////////////////////

for (iteration=0; iteration < DIRK_P700_number_of_tested_PARS; iteration++){
  
  //reset all important iteration-specific variables!
  var DIRK_P700_relaxation_data_vs_time=[];
  
  // now process each subtrace
  
  var DIRK_P700_trace_full = DIRK_P700_traces[iteration].data_raw; // get the ith trace that is measured
  
  var DIRK_P700_trace = DIRK_P700_trace_full.slice(DIRK_P700_trace_beginning, DIRK_P700_total_number_of_measurement_pulses); //ignore the first datapoints (often contain extreme values that are not interesting for us)
  
  ////////////////////////////////////
  // AVERAGE over all the subtraces // to get one single pulse snippet. 1) sum all subtraces, 2) divide by number_of_P700_subtraces
  ////////////////////////////////////
  
  var DIRK_P700_averaged_trace = DIRK_P700_trace.slice(0, DIRK_P700_length_of_subtrace); // the first subtrace must be added explicitly such that the following ones can be added
  
  for (var i = 1; i < DIRK_P700_number_of_subtraces; i++) {
      var next_subtrace = DIRK_P700_trace.slice(i * DIRK_P700_length_of_subtrace, (i+1) * DIRK_P700_length_of_subtrace);
      DIRK_P700_averaged_trace = TransformTrace('add', DIRK_P700_averaged_trace, next_subtrace);
  }	
  DIRK_P700_averaged_trace = TransformTrace('divide', DIRK_P700_averaged_trace, DIRK_P700_number_of_subtraces); // divide by the number of subtraces to get the average
  
  ///////////////////////////
  //remove spike artifacts //
  ///////////////////////////
  
  DIRK_P700_averaged_trace[DIRK_P700_length_of_baseline] = MathMEAN([DIRK_P700_averaged_trace[DIRK_P700_length_of_baseline-1], DIRK_P700_averaged_trace[DIRK_P700_length_of_baseline+1]]);
  DIRK_P700_averaged_trace[DIRK_P700_end_of_dark_period] = MathMEAN([DIRK_P700_averaged_trace[DIRK_P700_end_of_dark_period-1], DIRK_P700_averaged_trace[DIRK_P700_end_of_dark_period+1]]);
  
  DIRK_P700_averaged_trace[0] = MathMEAN([DIRK_P700_averaged_trace[1], DIRK_P700_averaged_trace[2]]);
  
  ////////////////////////////////////////////
  // Convert transmission to ABSORBANCE (A) //
  ////////////////////////////////////////////
  
  var DIRK_P700_baseline_value = MathMEAN(DIRK_P700_averaged_trace.slice(Math.floor(0.05*DIRK_P700_length_of_baseline), Math.floor(0.95*DIRK_P700_length_of_baseline)));
  // get the baseline trace from the baseline period in the averaged trace (avoid the first 5% and last 5% of the baseline because I am paranoid)
  
  // // in RIDES2.0 they use a linear regression over the kinetic part to find the baseline offset, but I see no goood reason as to why someone wouold do this??
  
  // infer the absorbance "abs" based on Beer's law (A = -log10(I/I0))
  var DIRK_P700_averaged_absorbtion = TransformTrace('abs', DIRK_P700_averaged_trace, DIRK_P700_baseline_value); // A
  
  // flipp the sign of every datapoint
  // var DIRK_P700_averaged_absorbtion = TransformTrace('subtract', DIRK_P700_averaged_absorbtion, TransformTrace('add', DIRK_P700_averaged_absorbtion, DIRK_P700_averaged_absorbtion)); // A
  
  ///////////////////////////////
  // PREPARE THE NONLINEAR FIT //
  ///////////////////////////////
  
  var DIRK_P700_dark_relaxation_window = DIRK_P700_averaged_absorbtion.slice(DIRK_P700_length_of_baseline-1, DIRK_P700_length_of_baseline + DIRK_P700_initial_dark_period_pulses_used_for_fitting);
  // when script was written, values were: (100-1, 120), we want 21 values in total
  
  // collect the timepoints and the values from the measurement into one object
  
  
  for (i = 0; i < DIRK_P700_initial_dark_period_pulses_used_for_fitting+1; i++) {
  	DIRK_P700_relaxation_data_vs_time.push([DIRK_P700_time_axis[i], DIRK_P700_dark_relaxation_window[i]]);
  }
  //output.fit_input = JSON.stringify(relaxation_data_vs_time); // for debugging and getting the window right
  
  ///////////////////////
  // NONLINEAR FITTING //
  ///////////////////////
  
  // now we fit the exponential decay to the function P700 = b + a * exp(-lambda * t), where we have a Baseline, Amplitude, and decay Constant
  // with lambda <- 1/c the nonlinear fit seems to perform MUCH MUCH better, so this is the formula in use below
  
  var DIRK_P700_fit = NonLinearRegression(DIRK_P700_relaxation_data_vs_time,{
    equation: decay,
    initial: [DIRK_P700_a_initial, DIRK_P700_b_initial, DIRK_P700_c_initial],
    iterations: 20000,
  });
  
  var DIRK_P700_a = DIRK_P700_fit.parameters[0].value;
  var DIRK_P700_b = DIRK_P700_fit.parameters[1].value;
  var DIRK_P700_c = DIRK_P700_fit.parameters[2].value;
  
  // output.a = a;
  // output.b = b; //this was helpful for debugging
  // output.c = c;
  
  // calculate the datapoints of the fit for optical comparison and sanity checks
  for (i = 0; i < DIRK_P700_initial_dark_period_pulses_used_for_fitting+1; i++) {
      var t = DIRK_P700_time_axis[i];
      var y = DIRK_P700_b + DIRK_P700_a * Math.exp(-t/DIRK_P700_c);
      DIRK_P700_fitted_curve.push(y);
  }
  
  //////////////////////
  // store everything //
  //////////////////////
  
  DIRK_P700_out_fit_R2.push(MathROUND(DIRK_P700_fit.r2      ,4));
  DIRK_P700_out_P700_amplitude.push(  MathROUND(DIRK_P700_a    ,4));      // transformed from absorbtion unit (AU) to milli absorbtion unit (mAU)
  DIRK_P700_out_k_PSI.push(        MathROUND(1 / DIRK_P700_c       ,4));      // unit: per millisecond (ms)
  DIRK_P700_out_P700_tau.push(        MathROUND(DIRK_P700_c           ,4));      // unit: millisecond
  DIRK_P700_out_v_PSI.push(        MathROUND((DIRK_P700_a*1000/ DIRK_P700_c) ,4));      // vH+ = P700_amplitude * gH+
  
  DIRK_P700_out_raw_transmission.push(DIRK_P700_trace);
  DIRK_P700_out_average_transmission.push(DIRK_P700_averaged_trace);
  DIRK_P700_out_averaged_absorbtion.push(DIRK_P700_averaged_absorbtion);
  DIRK_P700_out_relaxation_fit_data_and_fit.push([DIRK_P700_dark_relaxation_window, DIRK_P700_fitted_curve]);
  DIRK_P700_out_residuals_of_exponential_fit.push(TransformTrace("subtract", DIRK_P700_dark_relaxation_window, DIRK_P700_fitted_curve));
  DIRK_P700_out_dark_relaxation_window.push(DIRK_P700_dark_relaxation_window);
  DIRK_P700_out_fitted_curve_exp_decay = DIRK_P700_fitted_curve.concat(); //honestly, no clue why that works
  DIRK_P700_out_b.push(DIRK_P700_b); //sometimes interesting to look at
  
// LOOP END //
}           //
// LOOP END //


// FINALLY return all the outputs!

output.DIRK_P700_raw_transmission =       unwrap_if_single(DIRK_P700_out_raw_transmission);
output.DIRK_P700_average_transmission =   unwrap_if_single(DIRK_P700_out_average_transmission);
output.DIRK_P700_averaged_absorbtion =    unwrap_if_single(DIRK_P700_out_averaged_absorbtion);
output.DIRK_P700_relaxation_fit_data_and_fit = unwrap_if_single(DIRK_P700_out_relaxation_fit_data_and_fit); //only well-depicted in EDITOR and when only one PAR tested
// output.residuals_of_exponential_fit= unwrap_if_single(out_residuals_of_exponential_fit);
output.DIRK_P700_dark_relaxation_window =        unwrap_if_single(DIRK_P700_out_dark_relaxation_window);
output.DIRK_P700_fitted_curve_exp_decay=       unwrap_if_single(DIRK_P700_out_fitted_curve_exp_decay);
// output.b = unwrap_if_single(out_b)

output.DIRK_P700_fit_R2 =        unwrap_if_single(DIRK_P700_out_fit_R2);
output.DIRK_P700_amplitude_in_AU =   unwrap_if_single(DIRK_P700_out_P700_amplitude);  //milli absorbance units 
output.DIRK_P700_k_PSI_in_per_ms =         unwrap_if_single(DIRK_P700_out_k_PSI);        //millisecond
output.DIRK_P700_tau_in_ms =      unwrap_if_single(DIRK_P700_out_P700_tau);
output.DIRK_P700_v_PSI_in_AU_per_ms = unwrap_if_single(DIRK_P700_out_v_PSI);
output.DIRK_P700_time_between_pulses_in_ms = DIRK_P700_time_beween_pulses;



// ##########################################
// ###                                    ###
// ###        environmental vars          ###
// ###                                    ###
// ##########################################

var PAR = MathROUND(json.set[0].light_intensity,1);

var humidity1 = json.set[0].humidity;
var humidity2 = json.set[0].humidity2;
var humidity = MathROUND(MathMEAN([humidity1, humidity2]), 1);

var pressure1 = json.set[0].pressure;
var pressure2 = json.set[0].pressure2;
var pressure = MathROUND(MathMEAN([pressure1, pressure2]), 1);

var temperature1 = json.set[0].temperature;
var temperature2 = json.set[0].temperature2;
var env_temperature =  MathROUND(MathMEAN([temperature1, temperature2]), 1);

var leaf_temperature = MathROUND(json.set[0].contactless_temp, 1);

output.PAR = PAR;
output.env_temp = env_temperature;
output.leaf_temp = leaf_temperature;
output.rel_humidity = humidity;
output.pressure = pressure;

// Return Output Object (required)
return output;