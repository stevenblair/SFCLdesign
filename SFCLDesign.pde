/**
 * Visualisation of SFCL current-time characteristics
 *
 * Copyright (c) 2012 Steven Blair
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


final int ITERATIONS = 4096;
float i[] = new float[ITERATIONS];
float t[] = new float[ITERATIONS];
float[] majorTicksCurrent;
float[] majorTicksTime;

final boolean MODE_LINEAR = false;
final boolean MODE_LOGLOG = true;

final int MAX_WIDTH = 900;
final int MAX_HEIGHT = 600;
final int PLOT_START_X = 50;
final int PLOT_END_X = 850;
final int PLOT_START_Y = 550;
final int PLOT_HEIGHT = 530;

final int TICKS_PER_POWER_LINEAR_CURRENT = 2;
final int TICKS_PER_POWER_LINEAR_TIME = 4;

final color PlotColor = color(180, 33, 38, 255);  // red
final color IcColor = color(34, 177, 76, 200);	  // green

float sc_length;
float sc_radius;
float sc_area;
float Cv;
float kappa;
float Csc;
float Rsc;
float Jc_77K;
float E0;
float Ec;
float alpha_77K;
float beta;
float Tc;
float rho_Tc;
float Ta;
float T;

float I_PER_ITERATION, iMax, iMin;
float SCALE_CURRENT;
float SCALE_TIME;

public class Slider {
  float low;
  float high;

  public Slider(float low, float high) {
    this.low = low;
    this.high = high;
  }

  public float lowValue() {
    return low;
  }

  public float highValue() {
    return high;
  }
}

Slider currentRangeSlider = new Slider(0.0, 10000.0);
Slider timeRangeSlider = new Slider(0.0, 15.0);

void updateValues() {
  // superconductor dimensions and thermal properties
  sc_length = Processing.data.sc_length;            // superconductor length (m)
  sc_radius = Processing.data.sc_radius;           // superconductor radius (m)
  Cv = Processing.data.Cv;
  kappa = Processing.data.kappa;
  
  // superconductor characteristics
  Jc_77K = 1.5e7;               // critical current density at 77K (A/m^2)
  E0 = 0.1;                     // E-field for transition from superconducting state to flux-flow state (V/m)
  Ec = 1e-6 * 100;              // definition of E-field required for Jc (V/m)
  alpha_77K = 6.0;              // exponent value during superconducting state (dimensionless)
  beta = 3;                     // exponent value during flux-flow state (dimensionless)
  Tc = 95;                      // critical temperature (K)
  rho_Tc = Processing.data.rho_Tc;   // superconductor resistivity at Tc (ohm-m)
  Ta = Processing.data.Ta;
  T = Tc;

  majorTicksCurrent = getTicks(Processing.data.currentRange[0], Processing.data.currentRange[1], TICKS_PER_POWER_LINEAR_CURRENT);
  majorTicksTime = getTicks(Processing.data.timeRange[0], Processing.data.timeRange[1], TICKS_PER_POWER_LINEAR_TIME);
  
  iMax = max(majorTicksCurrent[majorTicksCurrent.length - 1], Processing.data.currentRange[1]);
  iMin = min(majorTicksCurrent[0], Processing.data.currentRange[0]);
  
  updateDerivedValues();
}

void updateDerivedValues() {
  sc_area = PI * sc_radius * sc_radius;                   // superconductor cross-sectional area (m^2)
  Csc = sc_length * sc_area * Cv;
  Rsc = 1.0 / (kappa * sc_length * 2 * sc_radius * PI);   // includes the total surface area (minus the "ends")

  I_PER_ITERATION = (iMax - iMin) / float(ITERATIONS);
}

int getCurrentPower(float maxCurrent) {
  int currentPower = 0;
  float current = maxCurrent;
  
  while (current > 9.999) {
    current /= 10.0;
    currentPower++;
  }
  
  return currentPower;
}

float[] getTicksOld(float maxCurrent, int ticksPerPower) {
  int currentPower = getCurrentPower(maxCurrent);
  int layers = floor(maxCurrent / pow(10, currentPower));
  
  if ((maxCurrent - layers * pow(10, currentPower)) / pow(10, currentPower) > 0.5) {
    layers++;
  }
  
  int totalTicks = layers * ticksPerPower;
  
  float[]  majorTicks = new float[totalTicks];
  
  for (int i = 0; i < layers; i++) {
    int tick = i * ticksPerPower;
    majorTicks[tick + 1] = (i + 1) * pow(10, currentPower);
    majorTicks[tick] = majorTicks[tick + 1] - (pow(10, currentPower) * 0.5);
  }
  
  return majorTicks;
}

float[] getTicksLinear(float low, float high, int ticksPerPower) {
  float x = high;
  int n = 0;
  float m = 1.0;
  
  if (x == 0.0) {
    return new float[]{0};
  }
  
  if (x < 1.0) {
    while (x < 1.0) {
      x = x * 10;
      n--;
    }
  }
  else if (x > 10.0) {
    while (x > 10.0) {
      x = x / 10;
      n++;
    }
  }
  else {
    n = 0;
  }
  m = x;
  
  int numberOfTicks = floor(x) * ticksPerPower + 1;  // add one to cater for low value
  int extra = 0;
  
  for (int i = 1; i < ticksPerPower; i++) {
    if (x - float(floor(x)) > float(i) / float(ticksPerPower)) {
      numberOfTicks++;
      extra++;
    }
  }
  
  float[] ticks = new float[numberOfTicks];
  
  if (Processing.data.linear == false) {
    ticks[0] = 0.01;
  }
  else {
    ticks[0] = 0.0;
  }
  
  for (int i = 1; i < numberOfTicks; i++) {
    ticks[i] = floor(m) * pow(10, n) * i / (numberOfTicks - extra - 1);
  }
  
  /*
  // trim ticks not visible (out of plot range) before creating array
  int start = 0;
  int end = numberOfTicks - 1;
  for (int i = 0; i < numberOfTicks; i++) {
    if (ticks[i] < low && i + 1 < numberOfTicks) {
      start = i + 1;
    }
    if (ticks[i] > high) {
      end = i;
      break;
    }
  }
  }*/
  
  /*if (high > 1000.0) {
    println(low + ", " + high);
    println(ticks);
    println();
  }*/
  
  //return subset(ticks, start, end - (start - 1));
  return ticks;
}

float[] getTicksLogLog(float low, float high) {
  int lowPower = floor(log10(low));
  int highPower = ceil(log10(high));
  
  float end = pow(10, highPower);
  float start = 0.01;
  if (low == 0.0) {
    start = 0.01;
  }
  else {
    start  = pow(10, lowPower);
  }
  
  int numberOfTicks = 1 + highPower - lowPower;
  
  if (numberOfTicks <= 1) {
    return new float[]{1.0, 10.0};
  }
  
  float[] ticks = new float[numberOfTicks];
  
  for (int i = 0; i < numberOfTicks; i++) {
    ticks[i] = pow(10, lowPower + i);
  }
  
  return ticks;
}

float[] getMinorLogLogTicks(float majorTicks[]) {
  int len = (majorTicks.length - 1) * 8;
  float[] minorTicks = new float[len];
  
  for (int i = 0; i < majorTicks.length - 1; i++) {
    for (int j = 0; j < 8; j++) {
      minorTicks[i * 8 + j] = majorTicks[i] * (j + 2);
    }
  }
  
  return minorTicks;
}

float[] getTicks(float low, float high, int ticksPerPower) {
  if (Processing.data.linear == false) {
    return getTicksLogLog(low, high);
  }
  else {
    return getTicksLinear(low, high, ticksPerPower);
  }
}
  
float plotPoint(int startOffset, float sign, float value, float minValue, float scaling) {
  if (Processing.data.linear == false) {
    return startOffset + sign * log10(value / minValue) * scaling;
  }
  else {
    return startOffset + sign * (value - minValue) * scaling;
  }
}

float plotPointY(int startOffset, float value, float minValue) {
  return plotPoint(startOffset, -1.0, value, minValue, SCALE_TIME);
}

float plotPointX(int startOffset, float value, float minValue) {
  return plotPoint(startOffset, 1.0, value, minValue, SCALE_CURRENT);
}
  
void plot(float[] xData, float[] yData, int c) {
  int start = 0;
  float maxTime = max(yData);

  // draw x-axis and y-axis
  stroke(50);
  strokeWeight(4);
  line(PLOT_START_X, PLOT_START_Y, PLOT_END_X, PLOT_START_Y);
  line(PLOT_START_X, PLOT_START_Y, PLOT_START_X, PLOT_START_Y - PLOT_HEIGHT);
  
  // draw major tick lines
  strokeWeight(1);
  colorMode(RGB);
  fill(255, 255, 255, 200);
  colorMode(HSB);
  textAlign(CENTER);
  
  if (Processing.data.linear == false) {
    SCALE_CURRENT = (PLOT_END_X - PLOT_START_X) / (majorTicksCurrent.length - 1);
    SCALE_TIME = (PLOT_HEIGHT) / (majorTicksTime.length - 1);
  }
  else {
    SCALE_CURRENT = float(PLOT_END_X - PLOT_START_X) / (Processing.data.currentRange[1] - Processing.data.currentRange[0]);
    SCALE_TIME = float(PLOT_HEIGHT) / Processing.data.timeRange[1] - Processing.data.timeRange[0];
  }
  
  for (int tick = 0; tick < majorTicksCurrent.length; tick++) {
    float xPoint = plotPointX(PLOT_START_X, majorTicksCurrent[tick], Processing.data.currentRange[0]);
    
    if (xPoint >= PLOT_START_X) {
      line(xPoint, PLOT_START_Y, xPoint, PLOT_START_Y - PLOT_HEIGHT);
      text(nf(majorTicksCurrent[tick] / 1000.0, 1, 2) + "kA", xPoint, PLOT_START_Y + 30);
    }
  }
  
  // draw minor tick lines
  strokeWeight(0.5);
  if (Processing.data.linear == false) {
    float[] minorTicksCurrent = getMinorLogLogTicks(majorTicksCurrent);
    
    for (int tick = 0; tick < minorTicksCurrent.length; tick++) {
      float xPoint = plotPointX(PLOT_START_X, minorTicksCurrent[tick], Processing.data.currentRange[0]);
      
      if (xPoint >= PLOT_START_X) {
        line(xPoint, PLOT_START_Y, xPoint, PLOT_START_Y - PLOT_HEIGHT);
      }
    }
    
    float[] minorTicksTime = getMinorLogLogTicks(majorTicksTime);
    
    for (int tick = 0; tick < minorTicksTime.length; tick++) {
      float yPoint = 0.0;
      yPoint = plotPointY(PLOT_START_Y, minorTicksTime[tick], Processing.data.timeRange[0]);
      
      if (yPoint <= PLOT_START_Y) {
        line(PLOT_START_X, yPoint, PLOT_END_X, yPoint);
      }
    }
  }
  
  strokeWeight(1);
  for (int tick = 0; tick < majorTicksTime.length; tick++) {
    float yPoint = plotPointY(PLOT_START_Y, majorTicksTime[tick], Processing.data.timeRange[0]);
    
    if (yPoint <= PLOT_START_Y) {
      line(PLOT_START_X, yPoint, PLOT_END_X, yPoint);

      float tickLabel;

      if (majorTicksTime[tick] >= 0.01) {
        tickLabel = nf(majorTicksTime[tick], 1, 2);
      }
      else {
        tickLabel = nf(majorTicksTime[tick], 1, 3);
      }

      text(tickLabel + "s", PLOT_START_X - 22, yPoint);
    }
  }
  
  // plot current-time curve
  strokeWeight(4);
  stroke(c, 255);
  noFill();
  beginShape(POLYGON);
  
  for (int x = 0; x < ITERATIONS; x++) {
    float xPoint = plotPointX(PLOT_START_X, xData[x], Processing.data.currentRange[0]);
    float yPoint = plotPointY(PLOT_START_Y, yData[x], Processing.data.timeRange[0]);
    
    if (yPoint >= PLOT_START_Y - PLOT_HEIGHT && yPoint <= PLOT_START_Y) {
      vertex(xPoint, yPoint);
    }
  }
  endShape();

  // find and show Ic (vertical asymptote)
  for (int x = 0; x < ITERATIONS; x++) {
    if (isNaN(yData[x]) || yData[x] < 0.0) {
      start = x;
      yData[x] = maxTime;
    }
    else {
      break;
    }
  }
  if (start >= 0 && start < ITERATIONS && xData[start] > iMin && xData[start] < iMax) {
    float xPoint = plotPointX(PLOT_START_X, xData[start], Processing.data.currentRange[0]);
    
    stroke(IcColor);
    fill(IcColor);
    textAlign(CENTER);
    strokeWeight(2);
    text("Ic = " + nf((xData[start]) / 1000.0, 1, 2) + "kA", xPoint, PLOT_START_Y + 15);
    line(xPoint, PLOT_START_Y, xPoint, PLOT_START_Y - PLOT_HEIGHT);
  }
}

void setup() {
  size(MAX_WIDTH, MAX_HEIGHT);
  background(0);

  font = createFont("SansSerif.plain", 13, true);
  textFont(font);
  textSize(12);
  textLeading(10);
  textAlign(CENTER);

  smooth();
  colorMode(HSB);
  frameRate(60);

  updateValues();

  reset();  // implemented in JavaScript
}

void draw() {
  if (Processing.data.change == true) {
    Processing.data.change = false;
    
    updateValues();

    background(0);
    textFont(font);
    //boolean foundE0 = false;
    
    for (int y = 0; y < ITERATIONS; y++) {
      i[y] = I_PER_ITERATION * y + iMin;
      float J = i[y] / sc_area;
      float k = i[y] * sc_length * E0 * pow((Ec / E0), (beta/alpha_77K)) * pow((J / Jc_77K), beta);

      float r = - (Csc*(-(Rsc*(Ta-Tc)*atan((-Ta-Tc+2 *Ta)/sqrt(-Ta*Ta+2 *Ta *Tc-Tc*Tc+4 *Tc*k *Rsc-308 *k* Rsc)))/sqrt(-Ta*Ta+2 *Ta *Tc-Tc*Tc+4 *Tc* k* Rsc-308 *k* Rsc)-1/2* Rsc *log(-Ta *(Ta+Tc)+Ta *Tc+(Tc-77) *k* Rsc+Ta*Ta)));
      t[y] = Csc*(-(Rsc*(Ta-Tc)*atan((-Ta-Tc+2 *T)/sqrt(-Ta*Ta+2 *Ta *Tc-Tc*Tc+4 *Tc*k *Rsc-308 *k* Rsc)))/sqrt(-Ta*Ta+2 *Ta *Tc-Tc*Tc+4 *Tc* k* Rsc-308 *k* Rsc)-1/2* Rsc *log(-T *(Ta+Tc)+Ta *Tc+(Tc-77) *k* Rsc+T*T))+r;
      
      /*if (isNaN(t[y])) {
        t[y] = -1.0;
      }*/

      /*float Jc = Jc_77K * ((Tc - Ta) / (Tc - 77));
      float alpha_alt = log10(E0 / Ec) / log10( pow((Jc_77K / Jc), (1 - (1/beta))) * pow((E0 / Ec), (1 / alpha_77K)) );
      float alpha1 = max(beta, alpha_alt);
      float E=Ec * pow((J / Jc), alpha1);
      
      if (foundE0 == false && E >= E0) {
        foundE0 = true;
        println(y + ", I0 = " + i[y]);
      }*/
    }

    plot(i, t, PlotColor);
  }
}

boolean isNaN(float f) {
  if (f > 0 || f < 0 || f == 0) {
    return false;
  }
  return true;
}

float log10(float x) {
  return (log(x) / log(10.0));
}

