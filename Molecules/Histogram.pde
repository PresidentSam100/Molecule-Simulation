class Histogram {
  int subdivisions;
  float l1, l2, u1, u2;
  float [] data;
  boolean init;
  int h, w;
  Histogram(int win, int hin, int sub) {
    h = hin;
    w = win;
    init = false;
    subdivisions = sub;
  }
  void update(float [] inData) {
    init = true;
    data = new float[inData.length];
    for (int i = 0; i < inData.length; i++) {
      data[i] = inData[i];
    }
  }
  void display() {
    if (!init) {return;}
    l1 = min(data);
    u1 = max(data);
    float subD = (u1-l1)/subdivisions;
    int [] heights = new int[subdivisions];
    stroke(0);
    line(0,0,0,h);
    line(0,h,w,h);
    text(nfc(l1,2),0,h+10);
    text(nfc(u1,2),w,h+10);
    for (int s = 0; s < subdivisions; s++) {
      int v = 0;
      for (int i = 0; i < data.length; i++) {
        if (data[i] > subD*s && data[i] < subD*(s+1)) {
          v++;
        }
      }
      heights[s] = v;
    }
    float rat = max(heights);
    for (int s = 0; s < subdivisions; s++) {
      fill(map(subD*(s+.5),.5*scaleV,1.5*scaleV,30,255),255,255);
      rect(1.0*w/subdivisions*s,h-heights[s]/rat*h,w*1.0/subdivisions,heights[s]/rat*h);
    }
    noFill();
  }
}
