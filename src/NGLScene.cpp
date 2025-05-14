#include "NGLScene.h"

#include <QMouseEvent>
#include <QGuiApplication>

#include "NGLScene.h"
#include <ngl/NGLInit.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/Transformation.h>
#include <ngl/Util.h>
#include <iostream>
#include <ngl/ShaderLib.h>

NGLScene::NGLScene()
{
  // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
  setTitle("Blank NGL");
}


NGLScene::~NGLScene()
{
  std::cout<<"Shutting down NGL, removing VAO's and Shaders\n";
}



void NGLScene::resizeGL(int _w , int _h)
{
  m_win.width  = static_cast<int>( _w * devicePixelRatio() );
  m_win.height = static_cast<int>( _h * devicePixelRatio() );
  m_project =  ngl::perspective(45.0f, float(m_win.width)/float(m_win.height), 0.01f, 200.0f);
}


void NGLScene::initializeGL()
{
  // we must call that first before any other GL commands to load and link the
  // gl commands from the lib, if that is not done program will crash
  ngl::NGLInit::initialize();
  glClearColor(0.7f, 0.7f, 0.7f, 1.0f);			   // Grey Background
  // enable depth testing for drawing
  glEnable(GL_DEPTH_TEST);
  // enable multisampling for smoother drawing
  glEnable(GL_MULTISAMPLE);
  ngl::VAOPrimitives::createSphere("sphere", 1.0f, 20);
  ngl::VAOPrimitives::createLineGrid("floor", 100, 100 ,50);
  m_emitter = std::make_unique<Emitter>(50000, 40000, 800, ngl::Vec3(0,0,0));
  ngl::ShaderLib::loadShader("ParticleShader", "shaders/ParticleVertex.glsl", "shaders/ParticleFragment.glsl");
  ngl::ShaderLib::use("ParticleShader");

  m_previousTime = std::chrono::steady_clock::now();

  m_text = std::make_unique<ngl::Text>("fonts/DejaVuSansMono.ttf",16);
  m_text -> setScreenSize(width(), height());
  m_text -> setColour(1,1,1);
  startTimer(10);

  m_view = ngl::lookAt({0,20,120}, {0,0,0}, {0,1,0});


}

void NGLScene::paintGL()
{
  // clear the screen and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0,0,m_win.width,m_win.height);
  ngl::ShaderLib::use("ParticleShader");
  auto rotX = ngl::Mat4::rotateX(m_win.spinXFace);
  auto rotY = ngl::Mat4::rotateY(m_win.spinYFace);
  auto mouseRotation = rotX * rotY;
  mouseRotation.m_m[3][0] = m_modelPos.m_x;
  mouseRotation.m_m[3][1] = m_modelPos.m_y;
  mouseRotation.m_m[3][2] = m_modelPos.m_z;

  glPointSize(5);
  ngl::ShaderLib::setUniform("MVP", m_project * m_view * mouseRotation);
  m_emitter -> render(m_win.width, m_win.height);

  //m_emitter -> debug();

  ngl::ShaderLib::use(ngl::nglColourShader);
  ngl::ShaderLib::setUniform("MVP", m_project * m_view * mouseRotation);
  ngl::ShaderLib::setUniform("Colour", 1.0f, 1.0f, 1.0f, 1.0f);
  ngl::VAOPrimitives::draw("floor");

  ngl::ShaderLib::use(ngl::nglTextShader);
  m_text -> renderText(10, 700, "Particle System");

}

// // Correct screenToWorld implementation
// ngl::Vec3 NGLScene::screenToWorld(int _x, int _y) const
// {
//   // Convert to NDC space [-1,1]
//   float x = (2.0f * _x) / width() - 1.0f;
//   float y = 1.0f - (2.0f * _y) / height();
//
//   // Create ray in clip space
//   ngl::Vec4 rayClip(x, y, -1.0f, 1.0f);
//
//   // Transform to eye space
//   ngl::Mat4 invProject = m_project;
//   invProject.inverse();
//   ngl::Vec4 rayEye = invProject * rayClip;
//   rayEye = ngl::Vec4(rayEye.m_x, rayEye.m_y, -1.0f, 0.0f);
//
//   // Transform to world space
//   ngl::Mat4 invView = m_view;
//   invView.inverse();
//   ngl::Vec4 rayWorld = invView * rayEye;
//   ngl::Vec3 rayDir(rayWorld.m_x, rayWorld.m_y, rayWorld.m_z);
//   rayDir.normalize();
//
//   // Camera position (assuming m_view is a lookAt matrix)
//   ngl::Vec3 eye(invView.m_30, invView.m_31, invView.m_32);
//
//   // Intersect with y=0 plane (ground)
//   if (rayDir.m_y != 0.0f)
//   {
//     float t = -eye.m_y / rayDir.m_y;
//     return eye + rayDir * t;
//   }
//
//   return ngl::Vec3(0.0f, 0.0f, 0.0f);
// }

//----------------------------------------------------------------------------------------------------------------------

void NGLScene::keyReleaseEvent(QKeyEvent* _event)
{
  m_keysPressed -= (Qt::Key)_event->key();
}


void NGLScene::keyPressEvent(QKeyEvent *_event)
{

  m_keysPressed += (Qt::Key)_event->key();
  // this method is called every time the main window recives a key event.
  // we then switch on the key value and set the camera in the GLWindow
  switch (_event->key())
  {

  case Qt::Key_M:
    m_emitter->toggleCursorInteraction();
    break;

  // escape key to quite
  case Qt::Key_Escape : QGuiApplication::exit(EXIT_SUCCESS); break;
  case Qt::Key_Space :

      m_win.spinXFace=0;
      m_win.spinYFace=0;
      m_modelPos.set(ngl::Vec3::zero());

  case Qt::Key_S:
    m_emitter->m_simulate = !m_emitter->m_simulate;
    break;


  break;
  case Qt::Key_A : m_animate ^= true; break;
  default : break;
  }
  // finally update the GLWindow and re-draw

    update();
}
void NGLScene::processKeys()
{
  float dx=0.0f;
  float dy=0.0f;
  float dz=0.0f;
  const float inc = 1.0f;
  for(auto key : m_keysPressed)
  {
    switch(key)
    {
      case Qt::Key_Left: dx -= inc; break;
      case Qt::Key_Right: dx += inc; break;
      case Qt::Key_Down: dz -= inc; break;
      case Qt::Key_Up: dz += inc; break;

    }
  }
  //m_emitter->move(dx, dy, dz);
}


void NGLScene::timerEvent(QTimerEvent *_event)
{
  auto now = std::chrono::steady_clock::now();
  auto delta = std::chrono::duration<float, std::chrono::seconds::period>(now - m_previousTime);
  m_previousTime = now;
  //std::cout <<"Time Delta" << delta.count() << "\n";
  if(m_animate)
  {
    processKeys();
    m_emitter -> update(delta.count());
  }
  update();
}
