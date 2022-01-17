///<music id="n9KjA471Y3k" loop="true"/>
///<texture id="https://media.istockphoto.com/photos/green-grass-field-picture-id182361617?b=1&k=20&m=182361617&s=170667a&w=0&h=2dn6dgL8aOXMC87GBJbQ7C4w1RzzLsTJEtE1zBL8bV8=">
uniform float u_time;
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform sampler2D u_texture0;
varying vec2 vUv;

#define MAX_STEPS 100 //Ray marching steps
#define MAX_DIST 200.
#define SURF_DIST .001 //Surface distance we have hit
#define LIGHT_COLOR vec3(0.782,0.91,0.965)

#define OBJECT_MASS 20.0
#define VELOCITY_X 0.0
#define VELOCITY_Y 10.0
#define VELOCITY_Z 0. 
#define PI 3.14159265359f
//---Volumetric clouds
const vec3 sunDir = normalize(vec3(-0.6f, 0.4f, 0.6f));
float hash( float n ) { return fract(sin(n)*43758.5453123); }
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
	
    float n = p.x + p.y*157.0 + 113.0*p.z;
    return 2.0f*mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                   mix( hash(n+157.0), hash(n+158.0),f.x),f.y),
               mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                   mix( hash(n+270.0), hash(n+271.0),f.x),f.y),f.z)-1.0f;
}
float fbm(in vec3 pos, int layers, float AM, float FM)
{
    float sum = 0.0f;
    float amplitude = 1.0f;
    for(int i=0;i<layers;++i)
    {
        sum += amplitude*noise(pos);
        amplitude *= AM;
        pos *= FM;
    }
    return sum;
}
float cloud(in vec3 p)
{
    return 0.01f*fbm(0.9f*vec3(0.2f, 0.2f, 0.3f)*(p+vec3(0.0f, 0.0f, .4f*u_time)), 7, 0.5f, 4.0f);
}
vec2 renderNoise(in vec3 ro, in vec3 rd)
{
    float tmin=10.0f;
    float tmax = 20.0f;
    float delta = 0.1f;
    float sum = 0.0f;
    float t = tmin;
    for(;t<tmax;t+=delta)
    {
        vec3 pos = ro + t*rd;
        //if(pos.y<-10.0f || pos.y>-1.0f || pos.x>5.0f || pos.x<-5.0f) return vec2(sum, -1.0f);
        
        float d = max(0.0f,cloud(pos));
        sum = sum*(1.0-d)+d;
        if(sum>0.99f) break;
    }
    return vec2(sum, t);
}
float shadeClouds(in vec3 ro, in vec3 rd)
{
    float sum =0.0f;
    float t = 0.0f;
    float delta = 0.1f;
    for(int i=0;i<5;++i)
    {
        vec3 pos = ro + rd*t;
        float d = max(0.0f,cloud(pos));
        sum = sum*(1.0-d)+d;
        if(sum>0.99f) break;
        t += delta;
    }
    return sum;
}
vec3 render(in vec3 ro, in vec3 rd)
{
    const vec3 sky = vec3(0.4, 0.6, 1.0);
    vec3 att = vec3(0.2f, 0.5f, 0.9f);
    vec2 ns = renderNoise(ro, rd);
    vec3 pos = ro+rd*ns.y;
    float shad = 1.0f;//0.9f*(1.0f-shadeClouds(pos+sunDir*0.1f, sunDir));
    float density = ns.x;
    float inv = (1.0f-density);
    
    float w = 1.8f*(0.5f*rd.y+0.5f);
    vec3 cl = shad*w*1.0f*mix(vec3(1.0f), inv*att, sqrt(density));
    if(density<0.1f) return mix(sky, cl, max(0.0f, density)*10.0f);
    //vec3 col = mix(sky, cl, 1.0f-exp(-0.0003*ns.y*ns.y) );
    //if(ns.y<0.0f) return sky;
	return cl; 
}
vec3 render(vec2 ndc, float aspectRatio)
{
	// camera origin
    vec3 o = vec3(0.0f, 0.0f, 0.0f);
	// camera horizontal field of view
    const float fov = 2.0f*PI / 3.0f;
    const float scaleX = tan(fov / 2.0f);
	// camera right vector
    vec3 right = vec3(1.0f, 0.0f, 0.0f)*scaleX;
	// camera forward vector
    vec3 forward = vec3(0.0f, 0.0f, 1.0f);
	// camera up vector
    vec3 up = vec3(0.0f, 1.0f, 0.0f)*scaleX*aspectRatio;
	// ray direction
    vec3 rd = normalize(forward + ndc.x*right + ndc.y*up);
    return render(o, rd);
}
//End of volumetric clouds

//Calculating physics position
vec3 PhysicsPos(float xVel, float upwardsVel,float zVel, float radius)
{
    float mass = OBJECT_MASS;
    vec3 force = vec3(0,-9.81,0);
    force += mass * vec3(0,-9.81,0);
    vec3 velocity = vec3(xVel,upwardsVel,zVel);
    velocity += force / mass * u_time;
    vec3 position = vec3(0,1,6);
    position += velocity * u_time;
   	
   	//Pseudo-collision
    if (position.y < radius)
    	position.y = radius;
   
    return position;
}
//Distance map of sphere
float SphereDist(vec3 p, vec3 pos, float radius)
{
	return length(p-pos)-radius;
}

float TorusDis(vec3 p, vec2 r)
{
	float x = length(p.xz)-r.x;
	return length(vec2(x, p.y))-r.y;
}

//Ray marching Distance from camera to object and plane
float GetDistance(vec3 p)
{
    //PHYSICS
    float radius = 1.03;
    vec3 Phys = PhysicsPos(VELOCITY_X, VELOCITY_Y, VELOCITY_Z, radius);
    
    //SPHERE RENDER
    vec4 sphere = vec4(Phys.xyz,radius); //Sphere position
    float sphereRay = length(p-sphere.xyz)-sphere.w;
    
    float dP = p.y;
    float d = min(sphereRay, dP);
    //d = min(d, SphereDist(p, vec3(Phys.x, Phys.y+1.3,Phys.z), 0.5*radius));
    
    //Torus render
    vec3 torusPosition = vec3(0,.6,6);
    d = min(d, TorusDis(p-torusPosition, vec2(1.5, 0.2)));
    
    for (float i = 1.; i < 5.; i++) {
      d = min(d, SphereDist(p, vec3(Phys.x, Phys.y+(i*0.55+.6),Phys.z), 0.05*i));
    }
    return d;
}

//Ray Marching
float RayMarch(vec3 ro, vec3 rd)
{
    float dO = 0.;
    for (int i = 0; i < MAX_STEPS; i++){
        vec3 p = ro + rd*dO;
        float dS = GetDistance(p);
        dO += dS;
        if (dO>MAX_DIST || dS < SURF_DIST) break;
    }
    return dO;
}
//PHONG LIGHTING
vec3 GetNormal(vec3 p)
{
    float d = GetDistance(p);
    vec2 e = vec2(.01,0);
    vec3 n = d - vec3(
    GetDistance(p-e.xyy),
    GetDistance(p-e.yxy),
    GetDistance(p-e.yyx));
    return normalize(n);
}
vec3 GetSpecularLighting(vec3 viewDir, vec3 reflectDir, float specularStrength)
{
    float spec1 = max(dot(viewDir, reflectDir), 0.0);
    float specTotal1 = 1.0;
    for (int i = 0; i < 64; i++)
        specTotal1 *= spec1;
    
    return specularStrength * specTotal1 * LIGHT_COLOR;
}
vec3 GetLight(vec3 p, vec3 ro)
{
    vec3 lightPos = vec3(0, 5, 6);
    lightPos.xz += vec2(sin(u_time), cos(u_time));
    vec3 l = normalize(lightPos - p);
    vec3 n = GetNormal(p);
    
    //Diffusion
    float dif = clamp(dot(n,l),0.,1.0);
    vec3 grassCol = texture(u_texture0,vUv).rgb;
    if (grassCol.g > 0.3 &&  grassCol.r > 0.3 && grassCol.b > 0.3)
    	grassCol.rgb = texture(u_texture0,vUv*0.32+0.3).rgb;
    vec3 diffusion = dif * grassCol;
    //Shadows
    float d = RayMarch(p+n*SURF_DIST*2.,l);
    if (d<length(lightPos-p)) diffusion *= 0.1;
    
    //Specular (BLINN-PHONG METHOD)
    float specularStrength = 0.5;
    vec3 viewDir = normalize(ro - p);
    vec3 reflectDir = reflect(-l, n);
    vec3 halfWayDir = normalize(l + viewDir);
    
    float spec = max(dot(viewDir, halfWayDir), 0.0);
    float specTotal = 1.0;
    for (int i = 0; i < 64; i++)
        specTotal *= spec;
       
    //Volumetric clouds
    
    vec3 volcloud = render(2.*gl_FragCoord.xy / u_resolution.xy-1., u_resolution.x/u_resolution.y);
    volcloud = sqrt(volcloud);
    vec3 BlinnPhong = specularStrength * specTotal * volcloud;
    
    vec3 specular1 = GetSpecularLighting(viewDir, reflectDir, specularStrength);
    
    return (diffusion + BlinnPhong + specular1);
}

void main()
{
    vec2 uv = (gl_FragCoord.xy-.5*u_resolution.xy)/u_resolution.y;
    
    vec3 col = vec3(0.0);
    
    //Camera model
    vec3 ro = vec3(0,1,0);
    //Ray direction (aka mouse rotation)
    //float xRotation = uv.x+0.03*u_mouse.x;
    vec3 rd = normalize(vec3(uv.x,uv.y,1.0));
    
    //Initializing Ray march
    float d = RayMarch(ro,rd);
    //Point to shade from lighting model
    vec3 p = ro + rd * d;
    vec3 diff = GetLight(p, ro);

    d /= 10.;
    col = diff;
    
    gl_FragColor = vec4(col,1.0);
}
