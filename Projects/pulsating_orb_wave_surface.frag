uniform vec2 u_resolution;
uniform float u_time;
const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float PRECISION = 0.001;

#define StarsNum 32.	//number of stars
#define StarsSize 0.035	//size of stars
#define StarsBright 2.0	//Bright of starss

struct Surface {
    float sd;
    vec3 col;
};

mat3 rotateY(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

Surface sdSphere(vec3 p, float r, vec3 offset, vec3 col)
{
  float d = length(p - offset) - r;
  d*=.3;
  vec3 gl = vec3(step(0.,-d));
  float glow = 0.01/d;
  
  return Surface(d, col);
}

Surface sdFloor(vec3 p, vec3 col) {
  float d = p.y + 1.;
  return Surface(d, col);
}

Surface minWithColor(Surface obj1, Surface obj2) {
  if (obj2.sd < obj1.sd) return obj2;
  return obj1;
}

Surface sdScene(vec3 p) {
  //remove "cos(u_time)" to remove pulsating effect
  Surface sphere = sdSphere(
  p, cos(u_time)+5.5+.5*sin(p.x*2.)*sin(p.z*2.)*sin(p.y*2.), vec3(0, 0, 0), vec3(.5,0,0)
  );
  
  return sphere;
}

Surface rayMarch(vec3 ro, vec3 rd, float start, float end) {
  float depth = start;
  Surface co;

  for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
    vec3 p = ro + depth * rd;
    co = sdScene(p);
    depth += co.sd;
    if (co.sd < PRECISION || depth > end) break;
  }
  
  co.sd = depth;
  
  return co;
}

vec3 calcNormal(in vec3 p) {
    vec2 e = vec2(1.0, -1.0) * 0.0005; // epsilon
    return normalize(
      e.xyy * sdScene(p + e.xyy).sd +
      e.yyx * sdScene(p + e.yyx).sd +
      e.yxy * sdScene(p + e.yxy).sd +
      e.xxx * sdScene(p + e.xxx).sd);
}
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
    return br*(smoothstep(.95, 1., (1. - sqrt(d))));
}


void main()
{
  vec2 uv = (gl_FragCoord.xy-.5*u_resolution.xy)/u_resolution.y;
  vec3 backgroundColor = vec3(stars(uv+u_time*0.03, StarsNum, StarsSize, StarsBright));

  vec3 col = vec3(0);
  vec3 ro = vec3(0, 0, 18);
  ro *= rotateY(u_time*.5);
  vec3 rd = normalize(vec3(uv, -1));
  rd*=rotateY(u_time*.5);
  
  Surface co = rayMarch(ro, rd, MIN_DIST, MAX_DIST);

  if (co.sd > MAX_DIST) {
    col = backgroundColor;
  } else {
    vec3 p = ro + rd * co.sd;
    vec3 normal = calcNormal(p);
    vec3 lightPosition = ro;
    vec3 lightDirection = normalize(lightPosition - p);

    // Calculate diffuse reflection by taking the dot product of 
    // the normal and the light direction.
    float dif = clamp(dot(normal, lightDirection), 0.5, 1.);

    // Multiply the diffuse reflection value by an orange color and add a bit
    // of the background color to the sphere to blend it more with the background.
    col = dif * co.col + vec3(0) * .2;
  }

  gl_FragColor = vec4(col, 1.0);
}
