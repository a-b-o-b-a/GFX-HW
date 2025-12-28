#version 300 es

precision highp float;
precision highp int;

struct Camera {
    vec3 pos;
    vec3 forward;
    vec3 right;
    vec3 up;
};
struct Ray
{
    vec3 origin;
    vec3 direction;
    
};

struct Plane {
    vec3 point;
    vec3 normal;
    vec3 color;
};

struct Sphere {
    vec3 center;
    float radius;
    vec3 color;
    int type; // 0: opaque, 1: reflective, 2: refractive
};

struct Light {
    vec3 position;
    vec3 direction;
    vec3 color;
    float shininess;
    float cutoff; // if > 0.0 then spotlight else directional light
};

struct HitInfo {
    vec3 rayOrigin;
    vec3 rayDir;
    float t;
    vec3 baseColor;
    int inside; // 1 if inside the sphere, 0 otherwise
    vec3 hitPoint;
    vec3 normal;
    bool isPlane;
    int type; // 0: diffuse, 1: reflective
};

const int TYPE_DIFFUSE = 0;
const int TYPE_REFLECTIVE = 1;
const int TYPE_REFRACTIVE = 2;

const int MAX_SPHERES = 16;
const int MAX_LIGHTS = 4;
const int MAX_DEPTH = 5;

in vec2 vUV;
out vec4 FragColor;

float EPSILON = 1e-3f;

uniform float uTime;
uniform ivec2 uResolution; // width and height of canvas


uniform Camera cam;
uniform Sphere uSpheres[MAX_SPHERES];
uniform int uNumSpheres;

uniform Light uLights[MAX_LIGHTS];
uniform int uNumLights;

uniform Plane uPlane;

vec3 checkerboardColor(in vec3 rgbColor, in vec3 hitPoint) {
    // Checkerboard pattern
    float scaleParameter = 2.0;
    float checkerboard = 0.0;
    if (hitPoint.x < 0.0) {
    checkerboard += floor((0.5 - hitPoint.x) / scaleParameter);
    }
    else {
    checkerboard += floor(hitPoint.x / scaleParameter);
    }
    if (hitPoint.z < 0.0) {
    checkerboard += floor((0.5 - hitPoint.z) / scaleParameter);
    }
    else {
    checkerboard += floor(hitPoint.z / scaleParameter);
    }
    checkerboard = (checkerboard * 0.5) - float(int(checkerboard * 0.5));
    checkerboard *= 2.0;
    if (checkerboard > 0.5) {
    return 0.5 * rgbColor;
    }
    return rgbColor;
}
vec3 at(float t, Ray ray) 
    {
        return ray.origin + t*ray.direction;
    }

bool intersectSphere(inout Ray ray, inout HitInfo hit, in Sphere sphere) 
    {
        vec3 rel = ray.origin-sphere.center;
        float a = dot(ray.direction, ray.direction);
        float b = 2.0f * dot(rel, ray.direction);
        float c = dot(rel, rel) - sphere.radius * sphere.radius;
        //find solution to quadratic equation
        float disc = b*b - 4.0*a*c;
        if (disc < 0.0) return false;
        //calculate both ts
        float t0 = (-b - sqrt(disc)) / (2.0*a);
        float t1 = (-b + sqrt(disc)) / (2.0*a);
        
        float t = -1.0f;
        // pick the closest positive t
        if (t0 > EPSILON) t = t0;
        else if (t1 > EPSILON) t = t1;
        else return false;
        //if actually hit, return the exact point to hit
        hit.rayOrigin = ray.origin;
        hit.rayDir = ray.direction;
        hit.type = sphere.type;
        hit.isPlane = false;
        hit.t = t;
        hit.hitPoint = at(t,ray);
        hit.normal = normalize(hit.hitPoint - sphere.center);
        hit.baseColor = sphere.color;
        return true;
    }
bool intersectPlane(inout Ray ray, inout HitInfo hit, in Plane plane) 
{
        float denom  = dot(plane.normal, ray.direction);
        if (abs(denom) < 1e-4) return false;

        float t = dot(plane.point - ray.origin, plane.normal) / denom;
        if (t <= 1e-3f) return false;  

    
        hit.rayOrigin = ray.origin;
        hit.rayDir = ray.direction;
        hit.isPlane = true;
        hit.t = t;
        hit.hitPoint = at(t, ray);
        //make sure the hit normal is pointing towards the light source, not at the same direction as the ray!
        if (dot(ray.direction, plane.normal) < 0.0)
        {
             hit.normal = normalize(plane.normal);
        }
        else
        {
            hit.normal = normalize(-(plane.normal));
        }
        hit.baseColor = plane.color;
        return true;
}

// replace void to your hit data type
// used for shadow calculation
bool intersectScene(inout Ray ray, float maxDist, inout HitInfo hit) 
{
    bool anyHit = false;

    hit.t = maxDist;
    for (int i=0;i<uNumSpheres;i++) {
        HitInfo temp;
        temp.t = 1e30f;
        if (intersectSphere(ray,temp,uSpheres[i]) && temp.t < maxDist && temp.t>EPSILON) {
            hit = temp;
            anyHit = true;
        }
    }
    
        HitInfo temp;
        temp.t = 1e30f;
        if (intersectPlane(ray,temp,uPlane) && temp.t < maxDist && temp.t>EPSILON) {
            hit = temp;
            anyHit = true;
        }
    

    return anyHit;
    

}

vec3 phong(inout HitInfo hit, vec3 viewDirection)
{
    //NO AMBIENT LIGHT??
    vec3 color =vec3(0.03,0.03,0.03);
    
    //directional lights
    for (int i=0;i<uNumLights;i++) {
        if(uLights[i].cutoff<0.0)
        {

        
        vec3 L = normalize(-uLights[i].direction);
        Ray shadowRay;
        HitInfo tempHit;
        shadowRay.origin = hit.hitPoint+ (hit.normal)*EPSILON;
        shadowRay.direction = L;
        if (!intersectScene(shadowRay, 1e30f,tempHit)) 
        {
        float diff = max(dot(hit.normal, L), 0.0f);

        vec3 R = reflect(-L, hit.normal);
        float spec = pow(max(dot(viewDirection, R), 0.0f), uLights[i].shininess);
        spec = min(spec, 1.0f);
        
        color +=  1.0 * 
            (hit.baseColor * diff +
            (0.7,0.7,0.7)* spec)
             * uLights[i].color;
        }
        }
        else
        //spot light

        {
             vec3 L = normalize(uLights[i].position - hit.hitPoint);
        
        Ray shadowRay;
        HitInfo tempHit;
        //offset to not hit the object itself
        shadowRay.origin = hit.hitPoint+ (hit.normal)*EPSILON;
        shadowRay.direction = L;

        float theta = dot(L, normalize(-uLights[i].direction));
        if ((theta > uLights[i].cutoff) && !intersectScene(shadowRay, length(uLights[i].position - hit.hitPoint) -EPSILON,tempHit)) 
        {
            float diff = max(dot(hit.normal, L), 0.0f);
            vec3 R = reflect(-L, hit.normal);
            float spec = pow(max(dot(viewDirection, R), 0.0f), uLights[i].shininess);

            //CLAMP SPECULAR! prevents insane levels of white close to camera on spots
            spec = min(spec, 1.0f);
            
            color += 1.0 * 
            ( hit.baseColor * diff +
             hit.baseColor * spec)
            * uLights[i].color;
        }
        }
    }

    if(hit.isPlane)
    {
        color = checkerboardColor(color, hit.hitPoint);
    }
    return clamp(color, 0.0f, 1.0f);
}

// replace int to your hit data type
/* calculates color based on hit data and uv coordinates */
vec3 calcColor(/*hit data type*/ int hitInfo) {

    return vec3(0);
}


vec3 raytrace(inout Ray ray) 
{
    vec3 result = vec3(0,0,0);
    bool lastRefracted= false;
    //iterate 5 times
    float reflectionPower = 1.0;
    for (int i=0;i<10;i++)
    {
    HitInfo closest;
    HitInfo tmp;
    closest.t = 1e30f;
    tmp.t = 1e30f;
    bool hitSomething = false;

    for (int i=0;i< uSpheres.length();i++)
    {
        if (intersectSphere(ray, tmp,uSpheres[i]) && tmp.t < closest.t)
            {
            closest = tmp;
            hitSomething = true;
            closest.isPlane = false;
        }
    }
        if (intersectPlane(ray,tmp,uPlane) && tmp.t < closest.t)
        {
            closest = tmp;
            hitSomething = true;
            closest.isPlane = true;
        }
    

    if (!hitSomething)
        {
            result = mix (result, vec3(0,0,0), reflectionPower);
            break;
        }
        // black background
        
    
    HitInfo hit = closest;
    vec3 viewDir = normalize(cam.pos - hit.hitPoint);
    vec3 localColor = phong(hit, viewDir);
    if(!lastRefracted) {
            result =mix (result, localColor, reflectionPower);
        } else {
            result = mix(result, localColor * reflectionPower, 0.7f);
        }
    if(hit.type==TYPE_DIFFUSE)
    {
        
        break;
    }
    else if(hit.type==TYPE_REFLECTIVE)
    {
        
        vec3 R = reflect(ray.direction, hit.normal);
        vec3 offset = dot(ray.direction, hit.normal) > 0.0f ? -hit.normal * EPSILON : hit.normal * EPSILON;
        ray.origin = hit.hitPoint + offset;
        ray.direction = normalize(R);
        reflectionPower *= 0.9f;
    }
    else if(hit.type==TYPE_REFRACTIVE)
    {
        
        float n1 = 1.0f;               
        float n2 = 1.5f;
        vec3 N = hit.normal;
        vec3 I = normalize(ray.direction);
        float cosI = dot(-I, N);
        float eta = n1 / n2;

        //flip normal if inside
        if (cosI < 0.0) {
            cosI = -cosI;
            N = -N;
            eta = n2 / n1;
        }

        float k = 1.0f - eta*eta*(1.0f - cosI*cosI);
        if (k >= 0.0f) 
        {
            vec3 T = eta * I + (eta * cosI - sqrt(k)) * N;
            ray.origin = hit.hitPoint + normalize(T) * EPSILON;
            ray.direction = normalize(T) ;
            
        } 
        else 
        {
            //internal reflection
            vec3 R = reflect(I, N);
            ray.origin = hit.hitPoint + N * EPSILON;
            ray.direction = normalize(R) ;
            
        }
        //some reflectance math??? for better result
        // Calculate Fresnel reflectance using original angle
        // float R0 = pow((n1 - n2) / (n1 + n2), 2.0);
        // float reflectance = R0 + (1.0 - R0) * pow(1.0 - abs(cosI), 5.0);
        result += localColor * 0.1f;
        // reflectionPower *= reflectance;
        lastRefracted=true;
         // attenuate for next bounce
        //float reflectance = 0.05 + 0.95 * pow(1.0 - cosI, 5.0); // Schlick approx
        //result += reflectionPower * localColor * (1.0 - reflectance);
    }
    //change the ray values to reflected ones, then calculate again 
    
    }
    
    
    return result;
    
}
/* scales UV coordinates based on resolution
 * uv given uv are [0, 1] range
 * returns new coordinates where y range [-1, 1] and x scales according to window resolution
 */
vec2 scaleUV(vec2 uv) {
    float scale = float(uResolution.x)/float(uResolution.y);
    float y = uv.y*2.0-1.0;
    float x = (uv.x*2.0-1.0)*scale;

    return vec2(x,y);
}

void main() {
    vec2 uv = scaleUV(vUV);
    vec3 rayDir = normalize(cam.forward + uv.x * cam.right + uv.y * cam.up);
    Ray ray ;
    ray.direction = rayDir;
    ray.origin = cam.pos;
    vec3 color = raytrace(ray);

    // test
    // remove when implementing raytracer
    //color = vec3(vUV.x, vUV.y, abs(sin(uTime)));
    FragColor = vec4(color, 1.0);
}

