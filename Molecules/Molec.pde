class Molec {
  PVector r, v, omega, ori;
  Matrix mi;
  float m, rad, l;
  ArrayList<PVector> centers;
  FloatList charges, masses;
  int cnum, mtype;
  Molec(PVector rin, PVector vin, float radin, float min, int type) {
    r = new PVector(rin.x,rin.y,rin.z);
    v = new PVector(vin.x,vin.y,vin.z);
    omega = new PVector();
    ori = PVector.random3D();
    mi = new Matrix(3,3);
    rad = radin;
    m = min;
    mtype = type;
    switch(type) {
      case 0: {
        cnum = 1;
        l = 0;
        centers = new ArrayList<PVector>(cnum);
        centers.add(new PVector());
        charges = new FloatList(cnum);
        charges.append(0.0);
        masses = new FloatList(cnum);
        masses.append(1);
        break;}
      case 1: {
        cnum = 2;
        l = .03;
        centers = new ArrayList<PVector>(cnum);
        centers.add(new PVector(-l,0,0));
        centers.add(new PVector(l,0,0));
        charges = new FloatList(cnum);
        charges.append(0.0);
        charges.append(0.0);
        masses = new FloatList(cnum);
        masses.append(1);
        masses.append(1);
        break;}
    }
    computeMomentOfInertiaTensor();
  }
  
  float getCharge(int i) {return charges.get(i);}
  PVector getC(int i) {return centers.get(i);}
  void computeMomentOfInertiaTensor() {
    //Ixx
    float hold = 0;
    for (int i = 0; i < cnum; i++) {
      hold += masses.get(i)*(sq(centers.get(i).y)+sq(centers.get(i).z));}
    mi.data[0][0] = hold;
    //Iyy
    hold = 0;
    for (int i = 0; i < cnum; i++) {
      hold += masses.get(i)*(sq(centers.get(i).x)+sq(centers.get(i).z));}
    mi.data[1][1] = hold;
    //Izz
    hold = 0;
    for (int i = 0; i < cnum; i++) {
      hold += masses.get(i)*(sq(centers.get(i).x)+sq(centers.get(i).y));}
    mi.data[2][2] = hold;
    //Ixy, Iyx
    hold = 0;
    for (int i = 0; i < cnum; i++) {
      hold -= masses.get(i)*centers.get(i).x*centers.get(i).y;}
    mi.data[0][1] = hold;
    mi.data[1][0] = hold;
    //Ixz, Izx
    hold = 0;
    for (int i = 0; i < cnum; i++) {
      hold -= masses.get(i)*centers.get(i).x*centers.get(i).z;}
    mi.data[0][2] = hold;
    mi.data[2][0] = hold;
    //Iyz, Izy
    hold = 0;
    for (int i = 0; i < cnum; i++) {
      hold -= masses.get(i)*centers.get(i).y*centers.get(i).z;}
    mi.data[1][2] = hold;
    mi.data[2][1] = hold;
  }
}
