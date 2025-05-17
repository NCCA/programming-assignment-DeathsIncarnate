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
                        this, &MainWindow::getValues);
            connect(m_ui->doubleSpinBox_2, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                        this, &MainWindow::getValues);
            connect(m_ui->doubleSpinBox_3, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                        this, &MainWindow::getValues);
            connect(m_ui->doubleSpinBox_4, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                        this, &MainWindow::getValues);
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
                m_boundingBox.m_width = static_cast<float>(m_ui->spinBox->value());
                m_boundingBox.m_depth = static_cast<float>(m_ui->spinBox_3->value());
                m_boundingBox.m_height = static_cast<float>(m_ui->spinBox_2->value());

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

void MainWindow::getValues(double value)
{
    try {
        std::cout << "Updating simulation parameters..." << std::endl;

        // Defensive checks
        if (!m_gl)
        {
            std::cerr << "Error: m_gl is null" << std::endl;
            return;
        }

        if (!m_ui)
        {
            std::cerr << "Error: m_ui is null" << std::endl;
            return;
        }

        // Get current values from all spinboxes with safety checks
        double yOffset = 0.0;
        double particleSpacing = 0.0;
        double pressureMultiplier = 0.0;
        double viscosityStrength = 0.0;

        try
        {
            yOffset = m_ui->doubleSpinBox->value();
            std::cout << "yOffset: " << yOffset << std::endl;

            particleSpacing = m_ui->doubleSpinBox_2->value();
            std::cout << "particleSpacing: " << particleSpacing << std::endl;

            pressureMultiplier = m_ui->doubleSpinBox_3->value();
            std::cout << "pressureMultiplier: " << pressureMultiplier << std::endl;

            viscosityStrength = m_ui->doubleSpinBox_4->value();
            std::cout << "viscosityStrength: " << viscosityStrength << std::endl;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Exception retrieving values: " << e.what() << std::endl;
            return;
        }

        // Call updateValues with verbose output
        std::cout << "Calling updateValues()..." << std::endl;
        m_gl->updateValues(pressureMultiplier, viscosityStrength, particleSpacing, yOffset);
        std::cout << "updateValues() completed successfully" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception in getValues: " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "Unknown exception in getValues" << std::endl;
    }
}