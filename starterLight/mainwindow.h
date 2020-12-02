#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <iostream>
#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/ArrayKernel.hh>
#include <OpenMesh/Core/Geometry/Vector11T.hh>

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
    // Fonction de Leo
    MyMesh::Point centre_gravite(MyMesh* _mesh);
    void boite_englobante(MyMesh* _mesh);
    std::map<int, int> valence(MyMesh* _mesh);
    // Fonction de Elias
    bool is_in_range(MyMesh::Scalar valueTest, MyMesh::Scalar a, MyMesh::Scalar marginOfError);
    void show_vf_list(MyMesh* _mesh);
    std::map<float, int> area_frequency(MyMesh* _mesh);
    std::map<MyMesh::Scalar, int> dihedral_angles(MyMesh* _mesh);
    void export_csv();
    void displayBox(MyMesh * _mesh, DisplayMode mode = DisplayMode::Normal);
    void delete_box(MyMesh * _mesh);
    // Fonctions de thomas
    float faceArea(MyMesh* _mesh, int faceID);
    float Ai(MyMesh* _mesh, int vertexID);
    float angleFF(MyMesh *_mesh, int faceID0, int faceID1, int vertID0, int vertID1);
    float angleEE(MyMesh* _mesh, int vertexID, int faceID);
    void H_Curv(MyMesh* _mesh);
    void K_Curv(MyMesh* _mesh);
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
    std::map<MyMesh::Scalar, int> ecart_ang(MyMesh* _mesh);
    bool is_hole(MyMesh* _mesh);


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

    void on_pushButton_clicked();

    void on_BoundingBox_clicked();

    void on_spinBox_valueChanged(int arg1);

private:

    bool modevoisinage;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    int n_e_old = 0;
    int n_f_old = 0;
    int n_v_old = 0;
    bool showBox = false;
    bool box_created = false;
    QString fileName;
    MyMesh::Scalar seuil = 1;
    Ui::MainWindow *ui;
};

template <typename T>
class MyStats
{
private:
    std::vector<T> _distrib ;
public:
    MyStats () {} ;

    void push_back (T data)
    {
        _distrib.push_back(data) ;
    }

    T min ()
    {
        T tmp = _distrib.at(0) ;
        for (int i=1 ; i<_distrib.size(); i++)
        {
            if (_distrib.at(i) < tmp)
                tmp = _distrib.at(i) ;
        }
        return tmp ;
    }

    T max ()
    {
        T tmp = _distrib.at(0) ;
        for (int i=1 ; i<_distrib.size(); i++)
        {
            if (_distrib.at(i) > tmp)
                tmp = _distrib.at(i) ;
        }
        return tmp ;
    }

    T mean ()
    {
        T acc(0) ;
        std::cout << "acc : " << acc << std::endl;
        if (_distrib.size() > 0)
        {
            for(int i=0; i<_distrib.size(); i++)
            {
                acc += _distrib.at(i) ;
            }
            return acc/(_distrib.size()) ;
        }
        else
            return acc ;
    }

    T stdev ()
    {
        T m = mean() ;
        T acc(0), tmp ;
        if (_distrib.size() > 0)
        {
            for(int i=0; i<_distrib.size(); i++)
            {
                tmp = _distrib.at(i) - m ;
                acc += tmp * tmp ;
            }
            return sqrt(acc/(_distrib.size())) ;
        }
        else
            return acc ;
    }

    T stdev (T m)
    {
        T acc, tmp ;
        if (_distrib.size() > 0)
        {
            for(int i=0; i<_distrib.size(); i++)
            {
                tmp = _distrib.at(i) - m ;
                acc += tmp * tmp ;
            }
            return sqrt(acc/(_distrib.size())) ;
        }
        else
            return acc ;
    }

    std::vector<T> get_vector(){
        return _distrib;
    }
};


#endif // MAINWINDOW_H
