#include "NGLScene.h"
#include <iostream>
#include <ngl/NGLInit.h>
#include <QMainWindow>
#include <ngl/VAOPrimitives.h>
#include <ngl/Transformation.h>
#include <ngl/Util.h>
#include <ngl/ShaderLib.h>
#include <Physics.h>


//----------------------------------------------------------------------------------------------------------------------
NGLScene::NGLScene( QWidget *_parent ) : QOpenGLWidget( _parent )
{

  // set this widget to have the initial keyboard focus
  setFocus();
  // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
  resize(_parent->size());
}

//----------------------------------------------------------------------------------------------------------------------
void NGLScene::initializeGL()
{

  ngl::NGLInit::initialize();
  glClearColor(0.4f, 0.4f, 0.4f, 1.0f);			   // Grey Background
  // enable depth testing for drawing
  glEnable(GL_DEPTH_TEST);

  // enable multisampling for smoother drawing
  glEnable(GL_MULTISAMPLE);
  ngl::VAOPrimitives::createSphere("sphere", 1.0f, 20);
  ngl::VAOPrimitives::createLineGrid("floor", 100, 100 ,50);
  emit glInitialized();

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

//----------------------------------------------------------------------------------------------------------------------
//This virtual function is called whenever the widget has been resized.
// The new size is passed in width and height.
void NGLScene::resizeGL( int _w, int _h )
{
  m_win.width  = static_cast<int>( _w * devicePixelRatio() );
  m_win.height = static_cast<int>( _h * devicePixelRatio() );
  m_project =  ngl::perspective(45.0f, float(m_win.width)/float(m_win.height), 0.01f, 200.0f);
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
    {
      m_emitter->toggleCursorInteraction();
      break;
    }

    // escape key to quite
  case Qt::Key_Escape : {QGuiApplication::exit(EXIT_SUCCESS); break;}
  case Qt::Key_Space :
    {
      m_win.spinXFace=0;
      m_win.spinYFace=0;
      m_modelPos.set(ngl::Vec3::zero());
    }

  case Qt::Key_S:
    {
      m_emitter->m_simulate = !m_emitter->m_simulate;
      break;
    }


    break;
  case Qt::Key_A : {m_animate ^= true; break;}
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
    case Qt::Key_Left: {dx -= inc; break;}
    case Qt::Key_Right: {dx += inc; break;}
    case Qt::Key_Down: {dz -= inc; break;}
    case Qt::Key_Up: {dz += inc; break;}

    }
  }
  //m_emitter->move(dx, dy, dz);
}

void NGLScene::updateValues(double currentPressure, double currentViscosity,
                          double currentParticleSpacing, double yOffset)
{
  m_emitter->m_physics->m_pressureMultiplier = static_cast<float>(currentPressure);
  m_emitter->m_physics->m_viscosityStrength = static_cast<float>(currentViscosity);
  m_emitter->m_physics->m_particleSpacing = static_cast<float>(currentParticleSpacing);
  m_emitter->yOffset = static_cast<float>(yOffset);

}

void NGLScene::simulationCheck()
{
    m_emitter->m_simulate = true;
}

void NGLScene::initialize()
{
  m_emitter->m_physics->m_maxParticles = 50000.0f;
  m_emitter->initializeParticles(m_emitter->yOffset, m_emitter->m_physics->m_particleSpacing, 50000);
  // m_emitter->m_boundingBox.updateVAO();
  // m_emitter->m_boundingBox.renderBoundingBox();
  // update();
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


NGLScene::~NGLScene()
{
}


