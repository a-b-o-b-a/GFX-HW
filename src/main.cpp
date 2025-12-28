
#include <glm/glm.hpp>
#include <Camera.h>
#include <fstream>
#include <iostream>
#include <vector>

#include <stb/stb_image.h>
#include <stb/stb_image_write.h>

using namespace std;
using namespace glm;

const int calcDepth = 20;
constexpr float EPSILON = 1e-3f;
const vec3 bgcolor = vec3(0,0,0);
struct Ray
{
    vec3 origin;
    vec3 direction;
    vec3 at(float t) const
    {
        return origin + t*direction;
    }
};

struct Material 
{
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    float shine;
    float alpha;
    float reflective = 0.0;   
};

struct RayCamera
{
    vec3 pos;
    vec3 up;
    vec3 direction;
    //default fov
    float fov = 70;
    //aspect = width/height, in this case 1
    float aspect;
    Ray get_ray(float u, float v) const {
        
        vec3 w = glm::normalize(direction);
        vec3 u_cam = glm::normalize(glm::cross(w, up));
        vec3 v_cam = glm::cross(u_cam, w);

        float theta = glm::radians(fov);
        float h = tan(theta * 0.5f);
        float viewport_height = 2.0f * h;
        float viewport_width  = aspect * viewport_height;

        vec3 horizontal = viewport_width * u_cam;
        vec3 vertical   = viewport_height * v_cam;

        vec3 lower_left = pos + w - horizontal * 0.5f - vertical * 0.5f;

        vec3 pixel = lower_left + u * horizontal + v * vertical;

        return {
            pos,
            normalize(pixel - pos)
        };
    }
};

struct GlobalLight 
{
    vec3 color;
    float intensity;
};

struct Light{
    virtual ~Light() {}
    vec3 color;
    float intensity;
    vec3 direction;
};
struct DirectionalLight : public Light
{
    
};

struct SpotLight : public Light
{  
    vec3 position; 
    float cutoff;     
};

struct Hit {
    float t;
    vec3 point;
    vec3 normal;
    //pointer to material of hit object
    Material* material;

    bool isPlane = false;
};
struct SceneObject{
    virtual ~SceneObject(){}
    Material material;
};
struct Sphere : public SceneObject
{
    vec3 center;
    float radius;
    
    bool intersect (const Ray& ray, Hit& hit) const
    {
        vec3 rel = ray.origin-center;
        float a = dot(ray.direction, ray.direction);
        float b = 2.0f * dot(rel, ray.direction);
        float c = dot(rel, rel) - radius * radius;
        //find solution to quadratic equation
        float disc = b*b - 4*a*c;
        if (disc < 0) return false;
        //calculate both ts
        float t0 = (-b - sqrt(disc)) / (2*a);
        float t1 = (-b + sqrt(disc)) / (2*a);
        
        float t = -1.0f;
        // pick the closest positive t
        if (t0 > EPSILON) t = t0;
        else if (t1 > EPSILON) t = t1;
        else return false;
        //if actually hit, return the exact point to hit
        hit.t = t;
        hit.point = ray.at(t);
        hit.normal = normalize(hit.point - center);
        hit.material = (Material*)&material;
        return true;
    }
};

struct Plane: public SceneObject
{
    vec3 position;
    vec3 normal;
    
    bool intersect(const Ray& ray, Hit& hit) const {
        float denom  = dot(normal, ray.direction);
        if (fabs(denom) < 1e-4) return false;

        float t = dot(position - ray.origin, normal) / denom;
        if (t <= 1e-3f) return false;  

    

        hit.t = t;
        hit.point = ray.at(t);
        //make sure the hit normal is pointing towards the light source, not at the same direction as the ray!
        if (dot(ray.direction, normal) < 0)
        {
             hit.normal = normalize(normal);
        }
        else
        {
            hit.normal = normalize(-normal);
        }
        hit.material = (Material*)&material;
        return true;
    }
};

struct Scene
{
    vector<Sphere> spheres;
    vector<Plane> planes;
    GlobalLight ambient;
    vector<SpotLight> spots;
    vector<DirectionalLight> dirLights;
    //check for any intersection below maxDist (where maxDist is distance to light source)
    bool intersect(const Ray& ray, float maxDist, Hit& hit) const {
    bool anyHit = false;
    constexpr float MIN_T = 1e-4f;
    hit.t = maxDist;
    for (const auto& sphere : spheres) {
        Hit temp;
        temp.t = 1e30f;
        if (sphere.intersect(ray, temp) && temp.t < maxDist && temp.t>MIN_T) {
            hit = temp;
            anyHit = true;
        }
    }
    for (const auto& plane : planes) {
        Hit temp;
        temp.t = 1e30f;
        if (plane.intersect(ray, temp) && temp.t < maxDist && temp.t>MIN_T) {
            hit = temp;
            anyHit = true;
        }
    }
    return anyHit;
}
};


vec3 checkerboardColor(vec3 rgbColor, vec3 hitPoint)
{
    // Checkerboard pattern
    float scaleParameter = 0.5f;
    float checkerboard = 0;
    if (hitPoint.x < 0)
        checkerboard += floor((0.5 - hitPoint.x) / scaleParameter);
    else
    {
        checkerboard += floor(hitPoint.x / scaleParameter);
    }
    if (hitPoint.y < 0)
    {
        checkerboard += floor((0.5 - hitPoint.y) / scaleParameter);
    }
    else
    {
        checkerboard += floor(hitPoint.y / scaleParameter);
    }
    checkerboard = (checkerboard * 0.5) - int(checkerboard * 0.5);
    checkerboard *= 2;
    if (checkerboard > 0.5)
    {
        return 0.5f * rgbColor;
    }
    return rgbColor;
}



vec3 phong(
    const Hit& hit,
    const vec3& viewDirection,
    const Scene& scene)
{
    vec3 color =scene.ambient.intensity * scene.ambient.color * hit.material->ambient;
    
    //directional lights
    for (const auto& light : scene.dirLights) {

        vec3 L = normalize(-light.direction);
        Ray shadowRay;
        Hit tempHit;
        shadowRay.origin = hit.point+ (hit.normal)*EPSILON;
        shadowRay.direction = L;
        if (!scene.intersect(shadowRay, 1e30f,tempHit)) {
        float diff = glm::max(dot(hit.normal, L), 0.0f);

        vec3 R = glm::reflect(-L, hit.normal);
        float spec = pow(glm::max(dot(viewDirection, R), 0.0f), hit.material->shine);
        spec = glm::min(spec, 1.0f);
        
        color += light.intensity * 
            (hit.material->diffuse * diff +
            hit.material->specular * spec)
             * light.color;
        }
    }

    //spot lights
    for (const auto& light : scene.spots) {
        vec3 L = normalize(light.position - hit.point);
        
        Ray shadowRay;
        Hit tempHit;
        //offset to not hit the object itself
        shadowRay.origin = hit.point+ (hit.normal)*EPSILON;
        shadowRay.direction = L;

        float theta = dot(L, normalize(-light.direction));
        if ((theta > light.cutoff) && !scene.intersect(shadowRay, length(light.position - hit.point) -EPSILON,tempHit)) {
            float diff = glm::max(dot(hit.normal, L), 0.0f);
            vec3 R = reflect(-L, hit.normal);
            float spec = pow(glm::max(dot(viewDirection, R), 0.0f), hit.material->shine);
            //CLAMP SPECULAR! prevents insane levels of white close to camera on spots
            spec = glm::min(spec, 1.0f);
            
            color += light.intensity * 
            (hit.material->diffuse * diff +
            hit.material->specular * spec)
            * light.color;
        }
    }
    if(hit.isPlane)
    {
        color = checkerboardColor(color,hit.point);
    }
    return clamp(color, 0.0f, 1.0f);
}


vec3 raytrace(const Ray& ray, int depth, const Scene& scene) {
    Hit closest;
    Hit tmp;
    closest.t = 1e30f;
    tmp.t = 1e30f;
    bool hitSomething = false;

    for (const auto& s : scene.spheres)
        if (s.intersect(ray, tmp) && tmp.t < closest.t)
            {
            closest = tmp;
            hitSomething = true;
            closest.isPlane = false;
        };

    for (const auto& p : scene.planes)
        if (p.intersect(ray, tmp) && tmp.t < closest.t)
        {
            closest = tmp;
            hitSomething = true;
            closest.isPlane = true;
        };
    

    if (!hitSomething)
        // black background
        return bgcolor;
    
    Hit hit = closest;
    vec3 viewDir = normalize(-ray.direction);
    vec3 localColor = phong(hit, viewDir, scene);
    vec3 result = localColor;
    if (hit.material->reflective>0 && depth < calcDepth )
    {
        vec3 R = reflect(ray.direction, hit.normal);
        //supposed to make sure the normal is pointing at the right direction for planes
        vec3 offset = dot(ray.direction, hit.normal) > 0
              ? -hit.normal * EPSILON
              : hit.normal * EPSILON;

        // //reflection glitch fix??
        //  // Adaptive offset for reflections
        // float cosTheta = abs(dot(hit.normal, R));
        // float reflectionBias = 1e-3f / glm::max(cosTheta, 0.1f);
        Ray reflected {
            hit.point + offset, 
            normalize(R)
        };

        vec3 reflectedColor = raytrace(reflected, depth + 1, scene);

        result = mix(result, reflectedColor, glm::min(hit.material->reflective,1.0f));
    }
    if (1.0-hit.material->alpha > 0 && depth < calcDepth) {
        //implement Snell's law
        float n1 = 1.0f;               
        float n2 = 1.5f;
        vec3 N = hit.normal;
        vec3 I = normalize(ray.direction);
        float cosI = dot(-I, N);
        float eta = n1 / n2;

        //flip normal if inside
        if (cosI < 0) {
            cosI = -cosI;
            N = -N;
            eta = n2 / n1;
        }

        float k = 1.0f - eta*eta*(1.0f - cosI*cosI);
        if (k >= 0.0f) 
        {
            vec3 T = eta * I + (eta * cosI - sqrt(k)) * N;
            Ray refracted { hit.point - N * EPSILON, normalize(T) };
            vec3 refractedColor = raytrace(refracted, depth + 1, scene);
            
            result = mix(result, refractedColor, 1-hit.material->alpha);
        } 
        else 
        {
            //internal reflection
            vec3 R = reflect(I, N);
            Ray reflected { hit.point + N * EPSILON, normalize(R) };
            vec3 reflectedColor = raytrace(reflected, depth + 1, scene);
            result = mix(result, reflectedColor, 1.0 - hit.material->alpha);
        }
       
    }

    return result;
    
}


void write_png(const char* filename, int width, int height, const std::vector<vec3>& framebuffer)
{
    std::vector<uint8_t> pixels(width * height * 3);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            vec3 c = framebuffer[y * width + x];
            int fy = height - 1 - y;
            int idx = 3 * (fy * width + x);

            pixels[idx+0] = (uint8_t)(255 * clamp(c.r, 0.0f, 1.0f));
            pixels[idx+1] = (uint8_t)(255 * clamp(c.g, 0.0f, 1.0f));
            pixels[idx+2] = (uint8_t)(255 * clamp(c.b, 0.0f, 1.0f));
        }
    }

    stbi_write_png(filename, width, height, 3, pixels.data(), width * 3);
}

vec3 getPlanePoint(float a, float b, float c, float d)
{
    vec3 point;
    if (c != 0)
        point = vec3(0, 0, -d / c);
    else if (b != 0)
        point = vec3(0, -d / b, 0);
    else if (a != 0)
        point = vec3(-d / a, 0, 0);
    
    return point;
}

void ProcessFile(FILE* file, Scene& scene, RayCamera& camera)
{
    string line;
    char t;
    float v1,v2,v3,v4;
    float distToScreen;
    float screenHeight;
    vector<Light*> sceneLights;
    vector<SceneObject*> sceneObjects;
    vector<SpotLight*> spots;
    int spotLightIndex = 0;
    int objectIndex = 0;
    int lightIntensityIndex = 0;
    while (fscanf(file," %c %f %f %f %f",&t,&v1,&v2,&v3,&v4) == 5)
    {
        switch (t)
        {
            case 'e':
                //camera position coordinates
                {
                    camera.pos = vec3(v1,v2,v3);
                    distToScreen = v4;
                    break;
                }
            case 'u':
                {
                    camera.up = vec3(v1, v2, v3);
                    screenHeight = v4;
                    camera.aspect = 1.0;
                    float fov = degrees(2 * atan(screenHeight / 2 / distToScreen));
                    camera.fov = fov;
                    break;
                }

            case 'f':
                camera.direction = vec3(v1, v2, v3);
                break;
            case 'a':
                {
                    GlobalLight ambient;
                    ambient.color = vec3(v1,v2,v3);
                    ambient.intensity = 1.0;
                    scene.ambient = ambient;
                    break;
                }
                
            case 'd':
                {
                if (v4==0.0)
                {
                    //directional light
                    DirectionalLight* dir = new DirectionalLight();
                    dir->direction = vec3(v1,v2,v3);
                    sceneLights.push_back(dir);
                    break;
                }
                else
                {
                    //spot light
                    SpotLight* spot = new SpotLight();
                    spot->direction = vec3(v1,v2,v3);
                    sceneLights.push_back(spot);
                    spots.push_back(spot);
                    break;
                    
                    
                }
                }
            case 'p':
                {
                    spots.at(spotLightIndex)->position = vec3(v1,v2,v3);
                    spots.at(spotLightIndex)->cutoff= v4; 
                    spotLightIndex++;
                    break;
                }
            case 'i':
                {

                sceneLights.at(lightIntensityIndex)->color = vec3(v1,v2,v3);
                sceneLights.at(lightIntensityIndex)->intensity = 1.0;
                lightIntensityIndex++;
                break;
                }

            case 'o':
                //opaque object
                if (v4>0)
                {
                    //sphere
                    Sphere* sp = new Sphere();
                    sp->center =  vec3(v1,v2,v3);
                    sp->radius = v4;
                    sceneObjects.push_back(sp);
                    sp->material.alpha = 1.0;
                    break;
                }
                else
                {
                    //plane
                    Plane* p = new Plane();
                    p->normal = normalize(vec3(v1,v2,v3));
                    p->position = getPlanePoint(v1,v2,v3,v4);
                    sceneObjects.push_back(p);
                    p->material.alpha = 1.0;
                    break;
                }
            case 'r':
                //reflective object, not yet implemented
                if (v4>0)
                {
                    //sphere
                    Sphere* sp = new Sphere();
                    sp->center =  vec3(v1,v2,v3);
                    sp->radius = v4;
                    sceneObjects.push_back(sp);
                    sp->material.alpha = 1.0;
                    sp->material.reflective = 1.0;
                    break;
                }
                else
                {
                    //plane
                    Plane* p = new Plane();
                    p->normal = normalize(vec3(v1,v2,v3));
                    p->position = getPlanePoint(v1,v2,v3,v4);
                    sceneObjects.push_back(p);
                    p->material.alpha = 1.0;
                    p->material.reflective = 1.0;
                    break;
                }
            case 't':
                //transparent object
                {

                
                if (v4>0)
                {
                    //sphere
                    Sphere* sp = new Sphere();
                    sp->center =  vec3(v1,v2,v3);
                    sp->radius = v4;
                    sceneObjects.push_back(sp);
                    sp->material.alpha = 0.3;
                    sp->material.reflective = 1.0;
                    break;
                }
                else
                {
                    //plane
                    Plane* p = new Plane();
                    p->normal = normalize(vec3(v1,v2,v3));
                    p->position = getPlanePoint(v1,v2,v3,v4);
                    sceneObjects.push_back(p);
                    p->material.alpha = 0.3;
                    p->material.reflective = 1.0;
                    break;
                }
                }
            case 'c':
                {
                
                sceneObjects.at(objectIndex)->material.ambient = vec3(v1,v2,v3);
                sceneObjects.at(objectIndex)->material.diffuse = vec3(v1,v2,v3);
                sceneObjects.at(objectIndex)->material.specular = vec3(0.7,0.7,0.7);
                sceneObjects.at(objectIndex)->material.shine = v4;
                objectIndex++;
                break;
                }
            default:
            
            break;
        }

    }
    for (SceneObject* obj: sceneObjects)
    {
        if (Sphere* s = dynamic_cast<Sphere*>(obj)) {
        scene.spheres.push_back(*s);
        }
        else if (Plane* p = dynamic_cast<Plane*>(obj)) {
        scene.planes.push_back(*p);
    }
    }
    for (Light* light: sceneLights)
    {
        if (DirectionalLight* l= dynamic_cast<DirectionalLight*>(light))
        {
            scene.dirLights.push_back(*l);
        }
        else if (SpotLight* s = dynamic_cast<SpotLight*>(light))
        {
            scene.spots.push_back(*s);
        }
    }
}
int main(int argc, char* argv[])
{
    
    int height = 1000;
    int width = 1000;
   
    for (int i=1;i<=6;i++)
    {
        Scene scene;
        RayCamera camera;
        vector<vec3> framebuffer(height*width);
        char filename[20];
        snprintf(filename, sizeof(filename), "./input/scene%d.txt", i);
        std::cout<<"started rendering scene: "<<i<<endl;
        FILE* file = fopen(filename, "r");
        ProcessFile(file,scene,camera);
        
        #pragma omp parallel for collapse(2)
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x)
            {

                float u = (x + 0.5f) / width;
                float v = (y + 0.5f) / height;

                Ray ray = camera.get_ray(u, v);
                vec3 color = raytrace(ray, 0, scene);

                framebuffer[y * width + x] = color;
            }
        }
        char outname[20];
        snprintf(outname, sizeof(outname), "./output/scene%d.png", i);
        std::cout<<"finished rendering scene: "<<i<<endl;
        write_png(outname,1000,1000,framebuffer);
        fclose(file);
        

    }
    for (int i=2;i<=5;i++)
    {
        Scene scene;
        RayCamera camera;
        vector<vec3> framebuffer(height*width);
        char filename[30];
        snprintf(filename, sizeof(filename), "./input/scene%d1.txt", i);
        std::cout<<"started rendering scene: "<<i<<"_1"<<endl;
        FILE* file = fopen(filename, "r");
        ProcessFile(file,scene,camera);
        
        #pragma omp parallel for collapse(2)
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x)
            {

                float u = (x + 0.5f) / width;
                float v = (y + 0.5f) / height;

                Ray ray = camera.get_ray(u, v);
                glm::vec3 color = raytrace(ray, 0, scene);

                framebuffer[y * width + x] = color;
            }
        }
        char outname[30];
        snprintf(outname, sizeof(outname), "./output/scene%d_1.png", i);
        std::cout<<"finished rendering scene: "<<i<<"_1"<<endl;
        write_png(outname,1000,1000,framebuffer);
        fclose(file);
        

    }
    
    return 0;
}