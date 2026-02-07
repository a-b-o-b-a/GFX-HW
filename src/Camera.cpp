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
            default:
                break;
        }
    }
}

void MouseButtonCallback(GLFWwindow* window, double currMouseX, double currMouseY)
{
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        std::cout << "MOUSE LEFT Click" << std::endl;
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

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        std::cout << "MOUSE LEFT Motion" << std::endl;
    }
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
    {
        std::cout << "MOUSE RIGHT Motion" << std::endl;
    }
}

void ScrollCallback(GLFWwindow* window, double scrollOffsetX, double scrollOffsetY)
{
    Camera* camera = (Camera*) glfwGetWindowUserPointer(window);
    if (!camera) {
        std::cout << "Warning: Camera wasn't set as the Window User Pointer! ScrollCallback is skipped" << std::endl;
        return;
    }

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