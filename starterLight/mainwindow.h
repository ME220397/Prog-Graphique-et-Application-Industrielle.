#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness; float value; Color faceShadingColor;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;


enum DisplayMode {Normal, TemperatureMap, ColorShading};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // les fonctions à compléter
    void get_carac(MyMesh* _mesh);
    bool onlyTriangles(MyMesh* _mesh);
    bool test_lonely_face(MyMesh* _mesh);
    bool test_lonely_vertex(MyMesh* _mesh);
    bool test_lonely_edge(MyMesh* _mesh);
    MyMesh::Normal face_norm(MyMesh* _mesh, MyMesh::FaceHandle face);
    MyMesh::Normal vertex_norm(MyMesh* _mesh, MyMesh::VertexHandle vertex);
    void print_vertices_norm(MyMesh* _mesh);
    void print_faces_norm(MyMesh* _mesh);
    void ecart_angulaire(MyMesh* _mesh);
    int nb_faces_isole(MyMesh* _mesh);
    int nb_points_isole(MyMesh* _mesh);
    int nb_aretes_isole(MyMesh* _mesh);

    void displayMesh(MyMesh *_mesh, DisplayMode mode = DisplayMode::Normal);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_angleArea_clicked();
    void on_pushButton_H_clicked();
    void on_pushButton_K_clicked();
    void on_pushButton_Affiche_Carac_clicked();
    void on_pushButton_normales_faces_clicked();
    void on_pushButton_normales_points_clicked();
    void on_pushButton_seuls_clicked();
    void on_pushButton_ecart_angulaire_clicked();

private:

    bool modevoisinage;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
