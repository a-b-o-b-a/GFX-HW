#include <stb/stb_image.h>
#include <stb/stb_image_write.h>
#include <glm/glm.hpp>
#include <iostream>
int main(void)
{
    std::string filepath = "Lenna.png";
    int width, height, comps;
    int req_comps = 4;
    unsigned char * buffer = stbi_load(filepath.c_str(), &width, &height, &comps, req_comps);
    
    //greyscale
    for( int i = 0;i<width*height;i++)
    {
        int sum = (buffer[i*comps+0] + buffer[i*comps+1] + buffer[i*comps+2]);
        //weighted pixels
        buffer[i*comps+0] = sum * 0.2989;
        buffer[i*comps+1] = sum * 0.5870;
        buffer[i*comps+2] = sum * 0.1140;
    }
    int result = stbi_write_png("new_lenna.png", width, height, req_comps, buffer, width * comps);
    std::cout << result << std::endl;
    return 0;
}