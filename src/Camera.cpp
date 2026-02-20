#include <Camera.h>
using namespace glm;
void Camera::SetOrthographic(float near, float far)
{
    m_Near = near;
    m_Far = far;

    // Rest Projection and View matrices
    m_Projection = glm::ortho(m_Left, m_Right, m_Bottom, m_Top, near, far);
    m_View = glm::lookAt(m_Position, m_Position + m_Orientation, m_Up);
}
void Camera::SetPerspective(float near, float far)
{
    m_Near = near;
    m_Far = far;

    // Rest Projection and View matrices
    m_Projection = glm::perspective(45.0f, 1.0f, near, far);
    m_View = glm::lookAt(m_Position, m_Position + m_Orientation, m_Up);
}
void Camera::SetPosition(glm::vec3 pos)
{
    m_Position = pos;
}
glm::vec3 Camera::GetPosition()
{
    return m_Position;
}
/////////////////////
// Input Callbacks //
/////////////////////

void KeyCallback(GLFWwindow* window, int key, int scanCode, int action, int mods)
{
    Camera* camera = (Camera*) glfwGetWindowUserPointer(window);
    if (!camera) {
        std::cout << "Warning: Camera wasn't set as the Window User Pointer! KeyCallback is skipped" << std::endl;
        return;
    }

    if (action == GLFW_PRESS || action == GLFW_REPEAT)
    {
        switch (key)
        {
            case GLFW_KEY_UP:
                std::cout << "UP Pressed" << std::endl;
                break;
            case GLFW_KEY_DOWN:
                std::cout << "DOWN Pressed" << std::endl;
                break;
            case GLFW_KEY_LEFT:
                std::cout << "LEFT Pressed" << std::endl;
                break;
            case GLFW_KEY_RIGHT:
                std::cout << "RIGHT Pressed" << std::endl;
                break;
            case GLFW_KEY_SPACE:
                std::cout << "setting clockwise: " << (camera->clockwise)<<std::endl;
                camera->clockwise = -1*(camera->clockwise);
                break;
            case GLFW_KEY_F:
                //rotate front wall
                std::cout << "rotating front wall " << camera->clockwise<<std::endl;
                for (Cube &cube : camera->cubes)
                {
                    //select front facing cubes
                    mat4 transform = mat4(1.0f);
                    transform = glm::translate(glm::mat4(1.0f), glm::vec3(cube.pos) ) * transform;
                    transform = cube.rot * transform;
                    vec4 worldCenter = transform * vec4(0.0f, 0.0f, 0.0f, 1.0f);
                    vec3 relPos = vec3(worldCenter);
                    //select those with z == 1
                    if (relPos.z >0.5)
                    {
                        mat4 newRot = glm::rotate(glm::mat4(1.0f), glm::radians(camera->rotAngle*(camera->clockwise)),vec3(0,0,1));
                        cube.rot = newRot*cube.rot;
                        std::cout << "rotating cube " << cube.pos.x<<cube.pos.y<<cube.pos.z<<std::endl;
                    }

                }
            break;
            case GLFW_KEY_B:
            //rotate back wall
            std::cout << "rotating back wall " << camera->clockwise<<std::endl;
            for (Cube &cube : camera->cubes)
            {
                //select front facing cubes
                mat4 transform = mat4(1.0f);
                transform = glm::translate(glm::mat4(1.0f), glm::vec3(cube.pos) ) * transform;
                transform = cube.rot * transform;
                vec4 worldCenter = transform * vec4(0.0f, 0.0f, 0.0f, 1.0f);
                vec3 relPos = vec3(worldCenter);
                
                if (relPos.z <-0.5)
                {
                    mat4 newRot = glm::rotate(glm::mat4(1.0f), glm::radians(camera->rotAngle*(camera->clockwise)),vec3(0,0,1));
                    cube.rot = newRot*cube.rot;
                    std::cout << "rotating cube " << cube.pos.x<<cube.pos.y<<cube.pos.z<<std::endl;
                }

            }
            break;
            case GLFW_KEY_R:
            
            std::cout << "rotating right wall " << camera->clockwise<<std::endl;
            for (Cube &cube : camera->cubes)
            {
                //select front facing cubes
                mat4 transform = mat4(1.0f);
                transform = glm::translate(glm::mat4(1.0f), glm::vec3(cube.pos) ) * transform;
                transform = cube.rot * transform;
                vec4 worldCenter = transform * vec4(0.0f, 0.0f, 0.0f, 1.0f);
                vec3 relPos = vec3(worldCenter);

                if (relPos.x >0.5)
                {
                    mat4 newRot = glm::rotate(glm::mat4(1.0f), glm::radians(camera->rotAngle*(camera->clockwise)),vec3(1,0,0));
                    cube.rot = newRot*cube.rot;
                    std::cout << "rotating cube " << cube.pos.x<<cube.pos.y<<cube.pos.z<<std::endl;
                }

            }
            break;
            case GLFW_KEY_L:

            std::cout << "rotating left wall " << camera->clockwise<<std::endl;
            for (Cube &cube : camera->cubes)
            {
                //select front facing cubes
                mat4 transform = mat4(1.0f);
                transform = glm::translate(glm::mat4(1.0f), glm::vec3(cube.pos) ) * transform;
                transform = cube.rot * transform;
                vec4 worldCenter = transform * vec4(0.0f, 0.0f, 0.0f, 1.0f);
                vec3 relPos = vec3(worldCenter);

                if (relPos.x <-0.5)
                {
                    mat4 newRot = glm::rotate(glm::mat4(1.0f), glm::radians(camera->rotAngle*(camera->clockwise)),vec3(1,0,0));
                    cube.rot =newRot*cube.rot;
                    std::cout << "rotating cube " << cube.pos.x<<cube.pos.y<<cube.pos.z<<std::endl;
                }

            }
            break;
            case GLFW_KEY_U:

            std::cout << "rotating up wall " << camera->clockwise<<std::endl;
            for (Cube &cube : camera->cubes)
            {
                //select front facing cubes
                mat4 transform = mat4(1.0f);
                transform = glm::translate(glm::mat4(1.0f), glm::vec3(cube.pos) ) * transform;
                transform = cube.rot * transform;
                vec4 worldCenter = transform * vec4(0.0f, 0.0f, 0.0f, 1.0f);
                vec3 relPos = vec3(worldCenter);

                if (relPos.y >0.5)
                {
                    mat4 newRot = glm::rotate(glm::mat4(1.0f), glm::radians(camera->rotAngle*(camera->clockwise)),vec3(0,1,0));
                    cube.rot = newRot*cube.rot;
                    std::cout << "rotating cube " << cube.pos.x<<cube.pos.y<<cube.pos.z<<std::endl;
                }

            }
            break;
            case GLFW_KEY_D:

            std::cout << "rotating down wall " << camera->clockwise<<std::endl;
            for (Cube &cube : camera->cubes)
            {
                //select front facing cubes
                mat4 transform = mat4(1.0f);
                transform = glm::translate(glm::mat4(1.0f), glm::vec3(cube.pos) ) * transform;
                transform = cube.rot * transform;
                vec4 worldCenter = transform * vec4(0.0f, 0.0f, 0.0f, 1.0f);
                vec3 relPos = vec3(worldCenter);

                if (relPos.y <-0.5)
                {
                    mat4 newRot = glm::rotate(glm::mat4(1.0f), glm::radians(camera->rotAngle*(camera->clockwise)),vec3(0,1,0));
                    cube.rot = newRot*cube.rot;
                    std::cout << "rotating cube " << cube.pos.x<<cube.pos.y<<cube.pos.z<<std::endl;
                }

            }
            break;

            case GLFW_KEY_A:
                if (camera->rotAngle < 180.0)
                {
                    camera->rotAngle = camera->rotAngle * 2;
                }
                std::cout << "setting rotation angle: " << camera->rotAngle << std::endl;
                break;
            case GLFW_KEY_Z:
            if (camera->rotAngle >90.0)
            {
                camera->rotAngle = camera->rotAngle / 2;
            }
            std::cout << "setting rotation angle: " << camera->rotAngle << std::endl;
            break;
            case GLFW_KEY_P:
            camera->colorPickMode = !(camera->colorPickMode);
            std::cout << "changing color pick mode to :  " << camera->colorPickMode << std::endl;
            break;
            default:
                break;
        }
    }
}

void MouseButtonCallback(GLFWwindow* window, double currMouseX, double currMouseY)
{
    Camera* camera = (Camera*) glfwGetWindowUserPointer(window);
    if (!camera) {
        std::cout << "Warning: Camera wasn't set as the Window User Pointer! KeyCallback is skipped" << std::endl;
        return;
    }
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        std::cout << "MOUSE LEFT Click" << std::endl;
        if(camera->colorPickMode)
        {
            std::cout << "picking color..." << std::endl;
            unsigned char color_picked[]{ 0, 0, 0, 0 };
            glReadPixels(currMouseX, camera->m_Height-currMouseY, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, color_picked);
            //find block
            for (Cube c : camera->cubes)
            {

                if(fabs(c.pickColor.x - float(color_picked[0])/255.0f)<0.01f && fabs(c.pickColor.y - float(color_picked[1])/255.0f)<0.01f &&fabs(c.pickColor.z - float(color_picked[2])/255.0f)<0.01f)
                {
                     std::cout << "cube picked:  " << c.pos.x<<' '<<c.pos.y<<' '<<c.pos.z<<std::endl;
                }
            }

        }
    }
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
    {
        std::cout << "MOUSE RIGHT Click" << std::endl;
    }
}

void CursorPosCallback(GLFWwindow* window, double currMouseX, double currMouseY)
{
    Camera* camera = (Camera*) glfwGetWindowUserPointer(window);
    if (!camera) {
        std::cout << "Warning: Camera wasn't set as the Window User Pointer! KeyCallback is skipped" << std::endl;
        return;
    }

    camera->m_NewMouseX = camera->m_OldMouseX - currMouseX;
    camera->m_NewMouseY = camera->m_OldMouseY - currMouseY;
    camera->m_OldMouseX = currMouseX;
    camera->m_OldMouseY = currMouseY;
    float rotationspeed = 0.01f;
    float movespeed = 0.01f;
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        if (camera->colorPickMode)
        {

            if (camera->pickedCubeId < 0)
            {
                std::cout << "picking color..." << std::endl;
                unsigned char color_picked[]{0, 0, 0, 0};
                glReadPixels(currMouseX, camera->m_Height - currMouseY, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, color_picked);
                // find block
                for (int i = 0; i < camera->cubes.size(); i++)
                {
                    Cube &c = camera->cubes[i];
                    if (fabs(c.pickColor.x - float(color_picked[0]) / 255.0f) < 0.01f && fabs(c.pickColor.y - float(color_picked[1]) / 255.0f) < 0.01f && fabs(c.pickColor.z - float(color_picked[2]) / 255.0f) < 0.01f)
                    {
                        std::cout << "cube picked:  " << c.pos.x << ' ' << c.pos.y << ' ' << c.pos.z << std::endl;
                        camera->pickedCubeId = i;
                    }
                }
            }
            else
            {
                // rotate picked cube
                float angleY = camera->m_NewMouseX * rotationspeed;
                mat4 rotY = rotate(mat4(1.0f), angleY, vec3(0.0f, 1.0f, 0.0f));

                vec3 right = normalize(cross(camera->m_Orientation, camera->m_Up));
                float angleX = camera->m_NewMouseY * rotationspeed;
                mat4 rotX = rotate(mat4(1.0f), angleX, right);
                camera->cubes[camera->pickedCubeId].rot = rotY * rotX * (camera->cubes[camera->pickedCubeId].rot);
            }
        }
        else
        {
            std::cout << "MOUSE LEFT Motion" << std::endl;


            float angleY = camera->m_NewMouseX * rotationspeed;
            mat4 rotY = rotate(mat4(1.0f), angleY, vec3(0.0f, 1.0f, 0.0f));


            vec3 right = normalize(cross(camera->m_Orientation, camera->m_Up));
            float angleX = camera->m_NewMouseY * rotationspeed;
            mat4 rotX = rotate(mat4(1.0f), angleX, right);

            vec4 newPos = rotY * rotX * vec4(camera->m_Position, 1.0f);
            camera->m_Position = vec3(newPos);

            camera->m_Orientation = normalize(-camera->m_Position);

            camera->m_View = lookAt(camera->m_Position, camera->m_Position + camera->m_Orientation, camera->m_Up);
        }
        
    }
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
    {
        std::cout << "MOUSE RIGHT Motion" << std::endl;

        vec3 right = normalize(cross(camera->m_Orientation, camera->m_Up));
        vec3 up = normalize(cross(right, camera->m_Orientation));

        camera->m_Position += right * float(camera->m_NewMouseX * movespeed);
        camera->m_Position += up * float(-camera->m_NewMouseY * movespeed);
        camera->m_View = lookAt(camera->m_Position, camera->m_Position + camera->m_Orientation, camera->m_Up);
    }
    else
    {
        camera->pickedCubeId = -1;
    }
}

void ScrollCallback(GLFWwindow* window, double scrollOffsetX, double scrollOffsetY)
{
    Camera* camera = (Camera*) glfwGetWindowUserPointer(window);
    if (!camera) {
        std::cout << "Warning: Camera wasn't set as the Window User Pointer! ScrollCallback is skipped" << std::endl;
        return;
    }
    float speed = 0.2f;
    float distance = length(camera->m_Position);
    float newDistance =  distance * (1.0-speed*scrollOffsetY);
    newDistance = clamp(newDistance,5.0f, 100.0f);

    camera->m_Position = normalize(camera->m_Position) * newDistance;
    camera->m_View = lookAt(camera->m_Position, camera->m_Position + camera->m_Orientation, camera->m_Up);
    std::cout << "SCROLL Motion" << std::endl;
}

void Camera::EnableInputs(GLFWwindow* window)
{
    // Set camera as the user pointer for the window
    glfwSetWindowUserPointer(window, this);

    // Handle key inputs
    glfwSetKeyCallback(window, (void(*)(GLFWwindow *, int, int, int, int)) KeyCallback);

    // Handle cursor buttons
    glfwSetMouseButtonCallback(window, (void(*)(GLFWwindow *, int, int, int)) MouseButtonCallback);

    // Handle cursor position and inputs on motion
    glfwSetCursorPosCallback(window , (void(*)(GLFWwindow *, double, double)) CursorPosCallback);

    // Handle scroll inputs
    glfwSetScrollCallback(window, (void(*)(GLFWwindow *, double, double)) ScrollCallback);
}