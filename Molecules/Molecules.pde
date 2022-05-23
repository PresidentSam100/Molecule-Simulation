import peasy.*;
PeasyCam cam;
ArrayList<Molec> gas;

HScrollbar hs1;
int lalign = 10;
int ualign = 18;
int spacer = 12;
int l = 220;
int r = 790;

float TEMPERATURE = 300;
int BOUNDARY = 0;
int TRACER = 0;
int MOLECULE = 0;
float scaleV;

//real constants
float kb = 1.3806e-23; //J/k
float an = 6.022e23;
//conversions to SI units
float tsr = 1e-12; //s; 100 fs, .1ps
float dr = 1e-9; //m ; 1 nm
float mr = 1.6605e-27; //kg; 1 dalton
float er = mr*dr*dr/tsr/tsr; //J; 1 janky joule
float kbr = kb/er;// jankyjoule/T

float massArgon = 39.948; //daltons, g/mol
float densityArgon = 1.922; //kg/m3
float radiusArgon = .188; //nm


float boxSize = 15; //nm
float density = densityArgon*1000/massArgon*an*pow(dr,3); // argon molecs/dr3 at 1ATM, 250K
int NUM = (int(density*pow(boxSize,3))/2)*2;
float nb = 1.333*3.1415*pow(radiusArgon,3)*NUM;

float effX;

int wid = 800;
int hei = 800;
int dep = 800;

int ts = 0;
float t = 0;

ArrayList<PVector> momentum;
ArrayList<PVector> prente;
ArrayList<PVector> tracer;
ArrayList<Float> idgRatio;
ArrayList<Float> vdwRatio;
Histogram speeds;
int storeLength = 200;
float ke = 0;
float pressureSum = 0;
int tracerIndex = NUM/2;
int tracerSize = 1000;

void setup() {
  hs1 = new HScrollbar(l,ualign,590,16,1);
  momentum = new ArrayList<PVector>();
  prente = new ArrayList<PVector>();
  tracer = new ArrayList<PVector>();
  idgRatio = new ArrayList<Float>();
  vdwRatio = new ArrayList<Float>();
  speeds = new Histogram(100,50,10);
  noFill();
  colorMode(HSB);
  size(800, 800, P3D);
  cam = new PeasyCam(this,400,400,400,400);
  gas = new ArrayList<Molec>(NUM);
  initializeGas();
  
}

void draw() {
  //visual environment
  background(0);
  
  drawP(wid,hei,dep);
  //data logging - makes new to hold for this cycle
  pressureSum = 0;
  prente.add(0,new PVector());
  momentum.add(0,new PVector());
  //collisions and movements
  for (int i = 0; i < NUM; i++){
    for (int j = i+1; j < NUM; j++) {
      isCollideBounceBalls(i,j);
    }
  }
  for (int i = 0; i < NUM; i++) {
    disp(i);
    move(i);
    checkBoundaries(i);
    //collect data
    momentum.set(0,PVector.add(momentum.get(0),gas.get(i).v)); //is overall momentum in dalton*nm/100fs
    prente.get(0).y += gas.get(i).v.magSq()*gas.get(i).m*.5; //in janky joules
    prente.get(0).z += gas.get(i).v.magSq()*gas.get(i).m/3.0/kbr; // in K
  }
  prente.get(0).x = pressureSum/6/boxSize/boxSize*NUM;
  if (momentum.size() == storeLength) {
    momentum.remove(storeLength-1);
    prente.remove(storeLength-1);
  }
  
  PVector prente1 = new PVector();
  PVector momentum1 = new PVector();
  for (int i =0; i < min(momentum.size(),storeLength); i++) {
    prente1.x += prente.get(i).x;
    prente1.y += prente.get(i).y;
    prente1.z += prente.get(i).z;
    momentum1 = PVector.add(momentum1,momentum.get(i));}
  prente1 = PVector.mult(prente1,1.0/momentum.size());
  momentum1 = PVector.mult(momentum1,1.0/momentum.size()/NUM);
  if (TRACER == 0) {
    tracer.add(0,new PVector(gas.get(tracerIndex).r.x,gas.get(tracerIndex).r.y,gas.get(tracerIndex).r.z));
    if (tracer.size() > tracerSize) {
      tracer.remove(tracerSize);
    }
    showTracer();}

  cam.beginHUD();
  noStroke();
  fill(255,0,255);
  rect(0,0,210,200);
  t = ts*tsr;
  fill(0);
  text("Time: \t" + nfc(t,13),lalign,ualign);
  text("Pressure: \t" + nfc(prente1.x, 6), lalign, ualign+1*spacer);
  text("Kinetic Energy: " + nfc(prente1.y,2),lalign,ualign+2*spacer);
  text("Temperature: " + nfc(prente1.z/2,2),lalign,ualign+3*spacer);
  text("Momentum: \t<" + nfc(momentum1.x,2)+", " + nfc(momentum1.y,2)+", "+nfc(momentum1.z,2)+">",lalign,ualign+4*spacer);
  float idgleft = prente1.x*pow(boxSize,3)/2;
  float vdwleft = prente1.x*(pow(boxSize,3)-nb)/2;
  float right = NUM*kbr*prente1.z/2;
  idgRatio.add(idgleft/right);
  vdwRatio.add(vdwleft/right);
  text("Left Side: \t" + nfc(idgleft,6),lalign,ualign+5*spacer);
  text("Right Side: \t" + nfc(right,6),lalign,ualign+6*spacer);
  float idgAvg = 0;
  float vdwAvg = 0;
  for (int i = 0; i < min(storeLength,idgRatio.size()); i++) {
    idgAvg+=idgRatio.get(i);
    vdwAvg+=vdwRatio.get(i);
  }
  idgAvg/=idgRatio.size();
  vdwAvg/=vdwRatio.size();
  text("IDG Ratio: \t" + nfc(idgAvg,4),lalign,ualign+7*spacer);
  text("VDW Ratio: \t" + nfc(vdwAvg,4),lalign,ualign+8*spacer);
  ts++;
  if (ts > storeLength-1) {
    idgRatio.remove(0);
    vdwRatio.remove(0);
  }
  if (ts % storeLength == 0 && ts > 0) {
    float[] indat = new float[NUM];
    for (int i = 0; i < NUM; i++) {
      indat[i] = gas.get(i).v.mag();
    }
    speeds.update(indat);
  }
  if (ts > storeLength+1) {
    pushMatrix();
      translate(lalign,ualign+9*spacer);
      speeds.display();
    popMatrix();}
  hs1.update();
  hs1.display();
  float low = .06;
  float high = .25;
  float helium = .14;
  float neon = .154;
  float argon = .188;
  float krypton = .202;
  float xenon = .216;
  float radon = .200;
  stroke(255,0,255);
  fill(255,0,255);
  line(map(helium,low,high,l,r),0,map(helium,low,high,l,r),40);
  line(map(neon,low,high,l,r),0,map(neon,low,high,l,r),40);
  line(map(argon,low,high,l,r),0,map(argon,low,high,l,r),40);
  line(map(krypton,low,high,l,r),0,map(krypton,low,high,l,r),40);
  line(map(xenon,low,high,l,r),0,map(xenon,low,high,l,r),40);
  line(map(radon,low,high,l,r),0,map(radon,low,high,l,r),40);
  pushMatrix();
    translate(map(helium,low,high,l,r),40);
    rotate(PI/2);
    text("Helium, 140 pm",0,0);
  popMatrix();
    pushMatrix();
    translate(map(neon,low,high,l,r),40);
    rotate(PI/2);
    text("Neon, 154 pm",0,0);
  popMatrix();
    pushMatrix();
    translate(map(argon,low,high,l,r),40);
    rotate(PI/2);
    text("Argon, 188 pm",0,0);
  popMatrix();
    pushMatrix();
    translate(5+map(krypton,low,high,l,r),40);
    rotate(PI/2);
    text("Krypton, 202 pm",0,0);
  popMatrix();
    pushMatrix();
    translate(map(xenon,low,high,l,r),40);
    rotate(PI/2);
    text("Xenon, 216 pm",0,0);
  popMatrix();
    pushMatrix();
    translate(map(radon,low,high,l,r),40);
    rotate(PI/2);
    text("Radon, 200 pm",0,0);
  popMatrix();
  noFill();
  
  if (map(hs1.getPos(),l,r,low,high) != radiusArgon) {
    radiusArgon = map(hs1.getPos(),l,r,low,high);
    for (int i = 0; i < NUM; i++) {
      gas.get(i).rad = radiusArgon/boxSize*wid;
    }
  }
  
  cam.endHUD();
  


}

void showTracer() {
  stroke(255,255,255);
  strokeWeight(4);
  beginShape();
  for (int k = 0; k < tracer.size(); k++) {
    vertex(tracer.get(k).x*wid/boxSize,tracer.get(k).y*wid/boxSize,tracer.get(k).z*wid/boxSize);}
  endShape();
  strokeWeight(1);
}

void initializeGas() {
  int unitCell = ceil(pow(NUM,.333));
  effX = boxSize;
  float gap = effX/(unitCell*1.); //1.07
  int n = 0;
  //building uniformly distributed velocity vectors
  scaleV = sqrt(3.*TEMPERATURE*kb/(massArgon*mr))/dr*tsr;
  println(scaleV);
  //float scaleO = sqrt(3.*TEMPERATURE*kb/(massArgon*mr))/dr*tsr;
  ArrayList<PVector>vs = new ArrayList<PVector>(NUM);
  ArrayList<PVector>os = new ArrayList<PVector>(NUM);
  for (int i = 0; i < NUM/2; i++) {
    float t = random(1)*6.283;
    float p = random(1)*3.142;
    float x1 = scaleV*cos(t)*sin(p);
    float y1 = scaleV*sin(t)*sin(p);
    float z1 = scaleV*cos(p);
    vs.add(new PVector(x1,y1,z1));
    vs.add(new PVector(-x1,-y1,-z1));
  }
  for (int xi = 0; xi < unitCell; xi++) {
    for (int yi = 0; yi < unitCell; yi++) {
      for (int zi = 0; zi < unitCell; zi++) {
        if (n < NUM) {
          float x = gap*(xi+.5);
          float y = gap*(yi+.5);
          float z = gap*(zi+.5);
          PVector rin = new PVector(x,y,z);
          PVector vin = new PVector(vs.get(n).x,vs.get(n).y,vs.get(n).z);
          gas.add(new Molec(rin, vin, radiusArgon/boxSize*wid, 1,MOLECULE));
          n++;     
  }}}}
  if (n < NUM-1) {
    NUM = n;
  }

}

void move(int i) {
  //TODO make sure obeys tennis racket theorem and conserves L
  gas.get(i).r.add(gas.get(i).v);
}

void isCollideBounceBalls(int i1, int i2) {
  PVector rd = PVector.sub(gas.get(i1).r,gas.get(i2).r);
  if (rd.mag() > 2*radiusArgon) {return;}
  float dll = (2*radiusArgon-rd.mag())/2;
  rd.setMag(dll);
  gas.get(i1).r.add(rd);
  gas.get(i2).r.add(PVector.mult(rd,-1));
  rd = PVector.sub(gas.get(i1).r,gas.get(i2).r);
  PVector vd = PVector.sub(gas.get(i1).v,gas.get(i2).v);
  float coeff11 = 2*gas.get(i2).m/(gas.get(i1).m+gas.get(i2).m);
  float coeff12 = PVector.dot(vd,rd)/rd.magSq();
  float coeff21 = 2*gas.get(i1).m/(gas.get(i1).m+gas.get(i2).m);
  float coeff22 = PVector.dot(PVector.mult(vd,-1),PVector.mult(rd,-1))/rd.magSq();
  PVector vec1 = new PVector(rd.x,rd.y,rd.z);
  PVector vec2 = new PVector(-rd.x,-rd.y,-rd.z);
  PVector vp1 = PVector.sub(gas.get(i1).v,PVector.mult(vec1,coeff11*coeff12));
  PVector vp2 = PVector.sub(gas.get(i2).v,PVector.mult(vec2,coeff21*coeff22));
  gas.get(i1).v = new PVector(vp1.x,vp1.y,vp1.z);
  gas.get(i2).v = new PVector(vp2.x,vp2.y,vp2.z);
}


void checkBoundaries(int index) {
    boolean crosser = false;
    switch (BOUNDARY) {
        case -1: { // no bounce
            break;}
        case 0: { //tiled
            if (gas.get(index).r.x > boxSize) {
                crosser = true;
                gas.get(index).r.x = gas.get(index).r.x - boxSize;
                pressureSum += gas.get(index).m*gas.get(index).v.x*2;}
            else if (gas.get(index).r.x < 0) {
                crosser = true;
                gas.get(index).r.x = gas.get(index).r.x + boxSize;
                pressureSum -= gas.get(index).m*gas.get(index).v.x*2;}
            if (gas.get(index).r.y > boxSize) {
                crosser = true;
                gas.get(index).r.y = gas.get(index).r.y - boxSize;
                pressureSum += gas.get(index).m*gas.get(index).v.y*2;}
            else if (gas.get(index).r.y < 0) {
                crosser = true;
                gas.get(index).r.y = gas.get(index).r.y + boxSize;
                pressureSum -= gas.get(index).m*gas.get(index).v.y*2;}
           if (gas.get(index).r.z > boxSize) {
                crosser = true;
                gas.get(index).r.z = gas.get(index).r.z - boxSize;
                pressureSum += gas.get(index).m*gas.get(index).v.z*2;}
            else if (gas.get(index).r.z < 0) {
                crosser = true;
                gas.get(index).r.z = gas.get(index).r.z + boxSize;
                pressureSum -= gas.get(index).m*gas.get(index).v.z*2;}        
            break;}
        }
        if (TRACER == 0) {
          if (index == tracerIndex && crosser) {
            tracer = new ArrayList<PVector>(tracerSize);
          }
        }
    }

void disp(int i) {
  pushMatrix();
  //println(scrnl*gas.get(i).r.x,scrnl*gas.get(i).r.y,scrnl*gas.get(i).r.z);
  translate(wid*gas.get(i).r.x/boxSize,wid*gas.get(i).r.y/boxSize,wid*gas.get(i).r.z/boxSize);
  stroke(cmap(gas.get(i).v.mag()),255,255);
  if (gas.get(i).centers.size() ==0) {
    sphere(gas.get(i).rad);
  }
  else {
    for (int j = 0; j < gas.get(i).centers.size(); j++) {
      pushMatrix();
      translate(wid*gas.get(i).centers.get(j).x*gas.get(i).ori.x/boxSize,
                wid*gas.get(i).centers.get(j).y*gas.get(i).ori.y/boxSize,
                wid*gas.get(i).centers.get(j).z*gas.get(i).ori.z/boxSize);
      sphere(gas.get(i).rad);
      popMatrix();
    }
  }
  popMatrix();
}

float cmap(float in) {
  return map(in,.5*scaleV,1.5*scaleV,30,255);}
void drawP(int x,int y,int z) {
  colorMode(HSB);
  stroke(255,0,255);
  beginShape();
    vertex(0,0,0);
    vertex(x,0,0);
    vertex(x,y,0);
    vertex(0,y,0);
    vertex(0,0,0);
    
    vertex(0,0,0);
    vertex(0,y,0);
    vertex(0,y,z);
    vertex(0,0,z);
    vertex(0,0,0);
    
    vertex(0,0,0);
    vertex(x,0,0);
    vertex(x,0,z);
    vertex(0,0,z);
    vertex(0,0,0);
    
    vertex(x,0,0);
    vertex(x,y,0);
    
    vertex(x,y,z);
    vertex(x,0,z);
    vertex(x,y,z);
    vertex(0,y,z);
    vertex(x,y,z);
    
    vertex(x,y,z);
    vertex(x,y,0);
    vertex(x,y,z);
    vertex(x,0,z);
    vertex(x,y,z);
    
    vertex(x,y,z);
    vertex(x,y,0);
    vertex(x,y,z);
    vertex(0,y,z);
    vertex(x,y,z);
  endShape();}    
