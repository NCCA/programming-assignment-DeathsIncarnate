#include "MainWindow.h"
#include "BoundingBox.h"
#include <iostream>

#include "ui_MainWindow.h"
#include "NGLScene.h"


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), m_ui(new Ui::MainWindow), m_gl(nullptr)
{
    try
    {
        m_ui->setupUi(this);

        // Initialize GL scene with proper error checking
        m_gl = new NGLScene(this);
        if (!m_gl)
            {
                std::cerr << "Error: Failed to create NGLScene" << std::endl;
                return;
            }

        // Add GL widget to layout
        m_ui->s_mainWindowGridLayout->addWidget(m_gl, 0, 0, 2, 1);


        // Connect all value spinboxes to update simulation parameters - use newer syntax with error checking
        try
        {
            connect(m_ui->doubleSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                        this, &MainWindow::updateSimulationParameters);
            connect(m_ui->doubleSpinBox_2, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                        this, &MainWindow::updateSimulationParameters);
            connect(m_ui->doubleSpinBox_3, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                        this, &MainWindow::updateSimulationParameters);
            connect(m_ui->doubleSpinBox_4, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                        this, &MainWindow::updateSimulationParameters);
            connect(m_ui->spinBox, QOverload<int>::of(&QSpinBox::valueChanged),
                    this, &MainWindow::updateSimulationParameters);
            connect(m_ui->spinBox_2, QOverload<int>::of(&QSpinBox::valueChanged),
                    this, &MainWindow::updateSimulationParameters);
            connect(m_ui->spinBox_3, QOverload<int>::of(&QSpinBox::valueChanged),
                    this, &MainWindow::updateSimulationParameters);
        }
        catch (const std::exception& e)
        {
                std::cerr << "Exception connecting signals: " << e.what() << std::endl;
        }

        // Connect bounding box update button
        connect(m_ui->pushButton_2, &QPushButton::pressed, this, [this]()
        {
            try
            {

                if (m_gl)
                {
                    m_gl->initialize();
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Exception in button handler: " << e.what() << std::endl;
            }
        });

        // Connect simulation toggle button
        try
        {
            connect(m_ui->pushButton, &QPushButton::clicked, m_gl, &NGLScene::simulationCheck);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Exception connecting simulation button: " << e.what() << std::endl;
        }

        std::cout << "MainWindow constructor completed successfully" << std::endl;
    }
    catch (const std::exception& e)
        {
            std::cerr << "Exception in MainWindow constructor: " << e.what() << std::endl;
        }
    catch (...)
        {
            std::cerr << "Unknown exception in MainWindow constructor" << std::endl;
        }
}

MainWindow::~MainWindow()
{
    delete m_ui;
}

void MainWindow::updateSimulationParameters()
{
    const double yOffset = m_ui->doubleSpinBox->value();
    const double particleSpacing = m_ui->doubleSpinBox_2->value();
    const double pressureMultiplier = m_ui->doubleSpinBox_3->value();
    const double viscosityStrength = m_ui->doubleSpinBox_4->value();
    const int width = m_ui->spinBox->value();
    const int height = m_ui->spinBox_2->value();
    const int depth = m_ui->spinBox_3->value();

    std::cout << "Updating parameters:\n"
              << "  yOffset: " << yOffset << "\n"
              << "  particleSpacing: " << particleSpacing << "\n"
              << "  pressureMultiplier: " << pressureMultiplier << "\n"
              << "  viscosityStrength: " << viscosityStrength << "\n"
              << "  Dimensions: " << width << "x" << height << "x" << depth << std::endl;

    m_gl->updateValues(pressureMultiplier, viscosityStrength, particleSpacing, yOffset, width, height, depth);
}