///<music id="mZvQ9ipTK_8" loop="true"/>
uniform vec2 u_resolution;
uniform float u_time;
const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float PRECISION = 0.001;
const float grid_intensity = 0.7;
#define PI 3.1415926
#define TAU 6.28318530718
#define MAX_ITER 5 //Water
#define StarsNum 32.	//number of stars
#define StarsSize 0.025	//size of stars
#define StarsBright 2.0	//Bright of stars
#define SPEED 15. //speed of camera

struct Surface {
    float sd; // signed distance value
    vec3 col; // color
};

//For instanced raymarching
Surface sdSphere(vec3 p, float r, vec3 offset, vec3 col)
{
  float d = length(mod(p,offset)-0.5*offset) - r; //modulo SDF instancing + adding '-0.5*c' for instancing
  return Surface(d, col);
}
//For regular object raymarching
Surface sdSphere1(vec3 p, float r, vec3 offset, vec3 col)
{
  float d = length(p - offset) - r;
  return Surface(d, col);
}

Surface sdFloor(vec3 p, vec3 col) {
  float d = p.y + 1.;
  return Surface(d, col);
}

Surface sdCappedCylinder(vec3 p, float h, float r, vec3 col){
  vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(h,r);
  float d1 = min(max(d.x,d.y),0.0) + length(max(d,0.0));
  return Surface(d1, col);
}

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h); }


Surface newRem(Surface s1, Surface s2, float strength) {
	float d1 = s1.sd;
	float d2 = s2.sd;
	float newD = opSmoothSubtraction(d1,d2, strength);
	return Surface(newD, s1.col);
}

Surface minWithColor(Surface obj1, Surface obj2) {
  if (obj2.sd < obj1.sd) return obj2; // The sd component of the struct holds the "signed distance" value
  return obj1;
}

Surface sdRoundBox( vec3 p, vec3 b, float r, vec3 col)
{
  vec3 q = abs(p) - b;
  float d = length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
  return Surface(d, col);
}

float pat(vec2 uv,float p,float q,float s,float glow)
{
    float z = cos(q * PI * uv.x) * cos(p * PI * uv.y) + cos(q * PI * uv.y) * cos(p * PI * uv.x);

    z += sin(u_time*4.0 + uv.x+uv.y * s)*.035;	// +wobble
    float dist=abs(z)*(1.0/glow);
    return dist;
}

//Water
vec3 water(vec2 uv) {
	float time = u_time * .5+23.0;
  vec2 p = mod(uv*TAU, TAU)-250.0;
	vec2 i = vec2(p);
	float c = 1.0;
	float inten = .005;
	for (int n = 0; n < MAX_ITER; n++) {
		float t = time * (1.0 - (3.5 / float(n+1)));
		i = p + vec2(cos(t - i.x) + sin(t + i.y), sin(t - i.y) + cos(t + i.x));
		c += 1.0/length(vec2(p.x / (sin(i.x+t)/inten),p.y / (cos(i.y+t)/inten)));
	}
	c /= float(MAX_ITER);
	c = 1.17-pow(c, 1.4);
	vec3 colour = vec3(pow(abs(c), 8.0));
    colour = clamp(colour + vec3(0.0, 0.35, 0.5), 0.0, 1.0);
    colour.b += 0.2;
    return colour;
}
//Sky
vec2 rand2(vec2 p)
{
    p = vec2(dot(p, vec2(12.9898,78.233)), dot(p, vec2(26.65125, 83.054543))); 
    return fract(sin(p) * 43758.5453);
}

float rand(vec2 p){
    return fract(sin(dot(p.xy ,vec2(54.90898,18.233))) * 4337.5453);
}

float stars(in vec2 x, float numCells, float size, float br)
{
    vec2 n = x * numCells;
    vec2 f = floor(n);

	float d = 1.0e10;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j){
            vec2 g = f + vec2(float(i), float(j));
			g = n - g - rand2(mod(g, numCells)) + rand(g);
            // Control size
            g *= 1. / (numCells * size);
			d = min(d, dot(g, g));
        }
    }

    return br * (smoothstep(.95, 1., (1. - sqrt(d))));
}
//City
float noise2d(vec2 p) {
	return fract(sin(dot(p ,vec2(17.9898,78.233))) * 456367.5453);
}
vec3 drawSky(){
	float resolution = max(u_resolution.y, u_resolution.y);
    vec2 coord = gl_FragCoord.xy / u_resolution;
    vec3 result = vec3(0.);
    
    //Sun
    vec2 pos = vec2(coord.x-0.5,coord.y-0.8);
    pos.y /= u_resolution.x/u_resolution.y;
    float dist = 1.0/length(pos);
    dist *= 0.1;
    dist = pow(dist, 0.9);
    vec3 col = dist * vec3(1.0, 0.5, 0.25);
    col = 1.0 - exp( -col );
    
    //City in back
    vec2 coord1 = 2.5*coord-0.75;
    float a = 0.;
    float size = 20.;
    for(float i = 0.; i < 1.; i+=1./size) {
		float t = size - i * size;
		coord1.x += i * 0.01;
		float px = floor((coord1.x + noise2d(vec2(i))) * t);
		float top = i * mix(i / size * sin(3.14 * 0.5 * 0.5), 0.7, 1.2);
		if(coord1.y + top < noise2d(vec2(px)))
			a = 1.0 * sin(3.14 * i);
	}
    vec3 city = vec3(a);
    
    result += stars(coord+0.002*u_time, StarsNum, StarsSize, StarsBright);
    result.r += 2.*-pos.y+0.35;
    result.b += 2.*-pos.y+0.35;
    return mix(city,mix(result,col,0.5)+0.1,0.65);
}

float grid(vec2 fragCoord, float space, float gridWidth)
{
    vec2 p  = fragCoord - vec2(.5);
    vec2 size = vec2(gridWidth - .5);
    
    vec2 a1 = mod(p - size, space);
    vec2 a2 = mod(p + size, space);
    vec2 a = a2 - a1;
       
    float g = min(a.x, a.y);
    return clamp(g, 0., 1.0);
}

Surface sdScene(vec3 p) {
  Surface sphereLeft = sdSphere1(p, 0., vec3(0.), vec3(0, .8, .8));
  Surface sphereRight = sdSphere1(p, 1., vec3(3,40,3), vec3(1, 0.58, 0.29));
  Surface co = minWithColor(sphereLeft, sphereRight);
  
  if (p.x < 6. && p.x > -6.){
  	co = minWithColor(co, sdCappedCylinder(mod(p,vec3(5.5,0.,15.5))-0.5*vec3(5.5,0.,15.5), .5, 5., vec3(1.)));
  	if (p.y < .0) {
  	vec3 pilOff = vec3(5.5,2.0,15.5);
  	co = minWithColor(co, sdRoundBox(mod(p,pilOff)-0.5*pilOff,
  				vec3(.5,0.025,.5), .25, vec3(1.)));
  	}
  	if (p.y < 10.){
  	
  	vec3 pilOff = vec3(5.5,10.0,15.5);
    vec3 pilOff1 = vec3(5.5,10.5,15.5);
  	Surface s1 = sdRoundBox(mod(p,pilOff)-0.5*pilOff,
  				vec3(.5,0.025,.5), .25, vec3(0.));
  	Surface s2 = sdRoundBox(mod(p,pilOff1)-0.5*pilOff1,
  				vec3(.5,0.025,.5), .25, vec3(1.));
  	Surface newSurf = newRem(s2,s1,.20);
  	co = minWithColor(co, newSurf);
  	}
  }
  vec3 roadCol = vec3(0.);
  vec3 checkerBoard = vec3(1. + 0.7*mod(floor(p.x) + floor(p.z), 2.0));
  
  float dF = pat(0.25*cos(p.zx), 5.0, 2.0, 35.0, 0.35);
  vec3 coolPattern =  vec3(0.925,0.45,1.25)*0.85/dF;
  
  //Checkerboard
  vec2 ass = p.xz;
  vec2 ass1 = p.xz /2.;
  vec3 gCol = coolPattern;
  gCol *= (1.-length(ass1-ass)/u_resolution.x*0.5);
  gCol *= grid(6.25*p.xz+0.5, 1., .2);

  roadCol = p.x > -2. && p.x < 2. ? gCol : water(0.5*vec2(p.x,p.z));
  vec3 floorColor = roadCol;
  co = minWithColor(co, sdFloor(p, floorColor));
  return co;
}

Surface rayMarch(vec3 ro, vec3 rd, float start, float end) {
  float depth = start;
  Surface co; // closest object

  for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
    vec3 p = ro + depth * rd;
    co = sdScene(p);
    depth += co.sd;
    if (co.sd < PRECISION || depth > end) break;
  }
  
  co.sd = depth;
  
  return co;
}
//Important for lighting calculations
vec3 calcNormal(in vec3 p) {
    vec2 e = vec2(1.0, -1.0) * 0.0005; // epsilon
    return normalize(
      e.xyy * sdScene(p + e.xyy).sd +
      e.yyx * sdScene(p + e.yyx).sd +
      e.yxy * sdScene(p + e.yxy).sd +
      e.xxx * sdScene(p + e.xxx).sd);
}

void main()
{
  vec2 uv = (gl_FragCoord.xy-.5*u_resolution.xy)/u_resolution.y;
  vec3 backgroundColor = drawSky();

  vec3 col = vec3(0);
  vec3 ro = vec3(0, 0, -SPEED*u_time); // ray origin that represents camera position
  vec3 rd = normalize(vec3(uv, -1)); // ray direction

  Surface co = rayMarch(ro, rd, MIN_DIST, MAX_DIST); // closest object

  if (co.sd > MAX_DIST) {
    col = backgroundColor; // ray didn't hit anything
  } else {
    vec3 p = ro + rd * co.sd; // point on sphere or floor we discovered from ray marching
    vec3 normal = calcNormal(p);
    vec3 lightPosition = vec3(ro.xy,ro.z-4.);
    vec3 lightDirection = normalize(lightPosition - p);

    float dif = clamp(dot(normal, lightDirection), 0.3, 1.);
    
    col = dif * co.col + vec3(0.05,0.2,0.84) * .2;
    col.r += 0.2;
  	col.b += 0.21;
  }
  gl_FragColor = vec4(col, 1.0);
}
