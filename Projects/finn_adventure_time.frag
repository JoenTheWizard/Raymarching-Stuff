///<music id="kKLzVv4hb8w" loop="true">
uniform vec2 u_resolution;
uniform float u_time;
uniform vec2 u_mouse;
const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float PRECISION = 0.001;
#define SHADES 5
//Background confs
#define TAU 6.28318530718
const vec3 BackColor	= vec3(0.0, 0.64, 0.78);
const vec3 CloudColor	= vec3(0.18,0.70,0.87);

struct Surface {
    float sd; // signed distance value
    vec3 col; // color
};

mat3 rotateX(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(1, 0, 0),
        vec3(0, c, -s),
        vec3(0, s, c)
    );
}

mat3 rotateY(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

mat3 rotateZ(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, -s, 0),
        vec3(s, c, 0),
        vec3(0, 0, 1)
    );
}


Surface sdSphere(vec3 p, float r, vec3 offset, vec3 col)
{
  float d = length(p - offset) - r;
  return Surface(d, col);
}

Surface sdTorus(vec3 p, vec2 t, vec3 col) {
	vec2 q = vec2(length(p.xz)-t.x,p.y);
	float d = length(q)-t.y;
	return Surface(d,col);
}

Surface sdBox(vec3 p, vec3 b, float r, vec3 col){
	vec3 q = abs(p) - b;
	float d = length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
	return Surface(d, col);
}

Surface sdCapsule( vec3 p, float h, float r, vec3 col){
  p.y -= clamp( p.y, 0.0, h );
  float d= length( p ) - r;
  return Surface(d, col);
}

Surface sdRoundBox(vec3 p, vec3 b, float r, vec3 col)
{
  vec3 q = abs(p) - b;
  float d =length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
  return Surface(d, col);
}

Surface sdFloor(vec3 p, vec3 col) {
  float d = p.y + 4.;
  return Surface(d, col);
}

Surface sdCappedCylinder(vec3 p, float h, float r, vec3 col)
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(h,r);
  float d1 = min(max(d.x,d.y),0.0) + length(max(d,0.0));
  return Surface(d1, col);
}

Surface minWithColor(Surface obj1, Surface obj2) {
  if (obj2.sd < obj1.sd) return obj2; // The sd component of the struct holds the "signed distance" value
  return obj1;
}

float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); 
}

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h); 
}

Surface newSub(Surface s1, Surface s2, float strength) {
	float d1 = s1.sd;
	float d2 = s2.sd;
	float newD = opSmoothUnion(d1,d2, strength);
	return Surface(newD, s1.col);
}

Surface newRem(Surface s1, Surface s2, float strength) {
	float d1 = s1.sd;
	float d2 = s2.sd;
	float newD = opSmoothSubtraction(d1,d2, strength);
	return Surface(newD, s1.col);
}

Surface sdEllipsoid(vec3 p, vec3 r, vec3 col)
{
  float k0 = length(p/r);
  float k1 = length(p/(r*r));
  float d = k0*(k0-1.0)/k1;
  return Surface(d, col);
}

Surface sdScene(vec3 p) {
  Surface sphereLeft = sdSphere(p, .07, vec3(-0.65, 0.72, -2.3), vec3(0.));
  Surface sphereRight = sdSphere(p, .07, vec3(.65, 0.72, -2.3), vec3(0.));
  Surface co = minWithColor(sphereLeft, sphereRight);
  
  //Torus
  //vec3 torusPosition = vec3(cos(u_time),1.5,-3.);
  //co = minWithColor(co, sdTorus(p-torusPosition, vec2(1.5,0.2), vec3(0.5,1.,1.)));
  
  //Ears
  vec3 ear1 = vec3(-1.4,1.5,-3.5);
  vec3 ear2 = vec3(1.4,1.5,-3.5);
  
  Surface earM = sdCapsule(p-ear1, .7,.4, vec3(1.));
  Surface earM1 = sdCapsule(p-ear2, .7,.4, vec3(1.));
  
  //Box
  vec3 boxPos = vec3(0.,.5,-4.);
  
  Surface boxM = sdBox(p-boxPos, vec3(1.,0.5,0.3), 1., vec3(1.));
  Surface newBox = newSub(boxM, earM, .5);
  Surface newBox1 = newSub(newBox, earM1, 0.5);
  co = minWithColor(co, newBox1);
  
  vec3 boxPos1 = vec3(0.,.5,-3.5);
  co = minWithColor(co, sdBox(p-boxPos1, vec3(.5,0.15,0.2), 1., vec3(1., 0.85, 0.76)));
  
  Surface mouth = sdEllipsoid(p-vec3(0.,0.35,-2.3), vec3(.35,0.2,0.1), vec3(0.));
  Surface mouthMesh = sdEllipsoid(p-vec3(0.,0.45,-2.3), vec3(.35,0.2,0.1), vec3(0.));
  
  Surface mouthModel = newRem(mouthMesh, mouth, 0.);
  co = minWithColor(co, mouthModel);
  //co = minWithColor(co, mouth);
  //co = minWithColor(co,mouthMesh);
  
  //Teeth
  Surface tooth1 = sdRoundBox(p-vec3(-0.19,.36,-2.2), vec3(.015, 0.007, 0.01), .03, vec3(1.));
  Surface tooth2 = sdRoundBox(p-vec3(-0.09,.36,-2.2), vec3(.015, 0.007, 0.01), .03, vec3(1.));
  Surface tooth3 = sdRoundBox(p-vec3(0.16,.36,-2.2), vec3(.015, 0.007, 0.01), .03, vec3(1.));
  co = minWithColor(co, tooth1);
  co = minWithColor(co, tooth2);
  co = minWithColor(co, tooth3);
  
  //Tongue
  Surface tongue = sdEllipsoid(p-vec3(-0.04,0.21,-2.3), vec3(.15,0.05,0.1), vec3(1.,0.,.0));
  co = minWithColor(co, tongue);
  
  //Torso
  vec3 torsoPos = vec3(0.,-2.7,-4.);
  Surface torso = sdBox(p-torsoPos, vec3(.85,1.3,0.3), 1., vec3(0., 0.6, 0.796));
  co = minWithColor(co, torso);
  
  //Backpack
  vec3 armBPos = vec3(-1.75,-1.25,-3.);
  Surface armB = sdTorus(((p-armBPos)*rotateX(1.5))*rotateZ(1.4), vec2(0.52,0.217), vec3(0.31,.651,.153));
  co = minWithColor(co, armB);
  
  vec3 armBPos1 = vec3(1.75,-1.25,-3.);
  Surface armB1 = sdTorus(((p-armBPos1)*rotateX(1.5))*rotateZ(-1.4), vec2(0.52,0.217), vec3(0.31,.651,.153));
  co = minWithColor(co, armB1);
  
  //Arms
  vec3 shirtCol = vec3(0., 0.6, 0.796);
  
  vec3 armPos = vec3(-2.1,-1.63,-2.8);
  vec3 armPos_0 = vec3(-2.1,-1.,-2.8);
  Surface arm1 = sdCappedCylinder(p-armPos, 0.35, .4, vec3(0., 0.6, 0.796));
  Surface arm1_0 = sdEllipsoid(p-armPos_0, vec3(0.37), shirtCol);
  Surface a_arm1 = newSub(arm1,arm1_0, 0.17);
  co = minWithColor(co, a_arm1);
  vec3 r_armPos = vec3(-2.1,-2.53,-2.8);
  Surface r_arm = sdCappedCylinder(p-r_armPos, 0.35, .5, vec3(1., 0.85, 0.76));
  co = minWithColor(co, r_arm);
  
  vec3 armPos1 = vec3(2.1,-1.63,-2.8);
  vec3 armPos1_0 = vec3(2.1,-1.,-2.8);
  Surface arm2 = sdCappedCylinder(p-armPos1, 0.35, .4, shirtCol);
  Surface arm2_0 = sdEllipsoid(p-armPos1_0, vec3(0.37), shirtCol);
  Surface a_arm = newSub(arm2,arm2_0, 0.17);
  co = minWithColor(co, a_arm);
  
  vec3 r_armPos1 = vec3(2.1,-2.53,-2.8);
  Surface r_arm1 = sdCappedCylinder(p-r_armPos1, 0.35, .5, vec3(1., 0.85, 0.76));
  co = minWithColor(co, r_arm1);
  
  vec3 floorColor = vec3(0.612,0.992, 0.478); //vec3(1. + 0.7*mod(floor(p.x) + floor(p.z), 2.0));
  
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

vec3 calcNormal(in vec3 p) {
    vec2 e = vec2(1.0, -1.0) * 0.0005; // epsilon
    return normalize(
      e.xyy * sdScene(p + e.xyy).sd +
      e.yyx * sdScene(p + e.yyx).sd +
      e.yxy * sdScene(p + e.yxy).sd +
      e.xxx * sdScene(p + e.xxx).sd);
}

//Background
float Func(float pX) {
	return 0.6*(0.5*sin(0.1*pX) + 0.5*sin(0.553*pX) + 0.7*sin(1.2*pX));
}
float FuncR(float pX) {
	return 0.5 + 0.25*(1.0 + sin(mod(40.0*pX, TAU)));
}
float Layer(vec2 pQ, float pT) {
	vec2 Qt = 3.5*pQ;
	pT *= 0.5;
	Qt.x += pT;

	float Xi = floor(Qt.x);
	float Xf = Qt.x - Xi -0.5;

	vec2 C;
	float Yi;
	float D = 1.0 - step(Qt.y,  Func(Qt.x));

	// Disk:
	Yi = Func(Xi + 0.5);
	C = vec2(Xf, Qt.y - Yi ); 
	D =  min(D, length(C) - FuncR(Xi+ pT/80.0));

	// Previous disk:
	Yi = Func(Xi+1.0 + 0.5);
	C = vec2(Xf-1.0, Qt.y - Yi ); 
	D =  min(D, length(C) - FuncR(Xi+1.0+ pT/80.0));

	// Next Disk:
	Yi = Func(Xi-1.0 + 0.5);
	C = vec2(Xf+1.0, Qt.y - Yi ); 
	D =  min(D, length(C) - FuncR(Xi-1.0+ pT/80.0));

	return min(1.0, D);
}
vec3 cloudBackground(vec2 UV) {
	vec3 Color = BackColor;
	for(float J=0.0; J<=1.0; J+=0.2)
	{
		// Cloud Layer: 
		float Lt =  u_time*(0.5  + 2.0*J)*(1.0 + 0.1*sin(226.0*J)) + 17.0*J;
		vec2 Lp = vec2(0.0, 0.3+1.5*( J - 0.5));
		float L = Layer(UV + Lp, Lt);

		// Blur and color:
		float Blur = 4.0*(0.5*abs(2.0 - 5.0*J))/(11.0 - 5.0*J);

		float V = mix( 0.0, 1.0, 1.0 - smoothstep( 0.0, 0.01 +0.2*Blur, L ) );
		vec3 Lc=  mix( CloudColor, vec3(1.0), J);

		Color =mix(Color, Lc,  V);
	}
	return Color;
}
//End of background

void main()
{
  vec2 uv = (gl_FragCoord.xy-.5*u_resolution.xy)/u_resolution.y;
  vec3 backgroundColor = cloudBackground(vec2(uv.x,-uv.y)*5.+1.2);
  //vec3(0.67, .84, .90);

  vec3 col = vec3(0);
  vec3 ro = vec3(0, 0, 3); // ray origin that represents camera position
  vec3 rd = normalize(vec3(uv, -1)); // ray direction
  //Mouse rotation
  //rd *= rotateX(u_mouse.y*0.005);
  //rd *= rotateY(u_mouse.x*0.005);
  
  Surface co = rayMarch(ro, rd, MIN_DIST, MAX_DIST); // closest object

  if (co.sd > MAX_DIST) {
    col = backgroundColor; // ray didn't hit anything
  } else {
    vec3 p = ro + rd * co.sd; // point on sphere or floor we discovered from ray marching
    vec3 normal = calcNormal(p);
    vec3 lightPosition = vec3(2, 2, 7);
    vec3 lightDirection = normalize(lightPosition - p);

    // Calculate diffuse reflection by taking the dot product of 
    // the normal and the light direction.
    float dif = clamp(dot(normal, lightDirection), 0.42, .73);

    // Multiply the diffuse reflection value by an orange color and add a bit
    // of the background color to the sphere to blend it more with the background.
    vec3 ambientLight = vec3(1.55);
    col = dif * co.col + ambientLight * .2;
  }

  // Output to screen
  gl_FragColor = vec4(col, 1.0);
  
  //Cel shading
  float brightness = (gl_FragColor.r + gl_FragColor.g + gl_FragColor.g) / 3.; 
  float shade = floor(brightness * float(SHADES));
  float brighnessOfShade = shade / float(SHADES);
  float factor = brightness / brighnessOfShade;
  if (!(co.sd > MAX_DIST)) {
  gl_FragColor.rgb /= vec3(factor);
  }
}
