#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <QDir>


/* **** début de la partie à compléter **** */
void MainWindow::translate_to_origin(MyMesh *_mesh){
    MyMesh::Point B;
    MyMesh::Point O;
    B *= 0;
    O *= 0;

    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++){
        B += _mesh->point(*v);
    }
    B /= _mesh->n_vertices();

    MyMesh::Point BO = O - B;

    for (MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++) {
        MyMesh::Point newPoint = _mesh->point(*v) + BO;
        _mesh->set_point(*v, newPoint);
    }

}
void MainWindow::get_carac(MyMesh* _mesh){
    int euler_formula = _mesh->n_vertices() - _mesh->n_edges() + _mesh->n_faces();
    qDebug() << "Nombre de sommets "<< _mesh->n_vertices();
    qDebug() << "Nombre de faces :" << _mesh->n_faces();
    qDebug() << "Nombre d'arêtes" << _mesh->n_edges();
    qDebug() << "Is triangle" << _mesh->is_triangles();
    qDebug() << "Euler-Poincaré = " << euler_formula;
}

void MainWindow::export_csv(){
    std::map<double, int> area_freq = area_frequency(&mesh);
    std::map<MyMesh::Scalar, int> dihedral_freq = dihedral_angles(&mesh);
    std::map<uint, int> valence_freq = valence(&mesh);
    std::map<MyMesh::Scalar, int> ecart_angulaire = ecart_ang(&mesh);
    /** Frequence des aires **/
    // chemin à changer selon votre systeme
    // Elias = "/Users/eliasmunoz/Documents/Git Projects/Prog-Graphique-et-Application-Industrielle./CSV/area.csv";
    QString path =  "./CSV/area.csv";
    //QString path = "/home/kammerlocher/prog_indus/Prog-Graphique-et-Application-Industrielle./CSV/area.csv";
    //QString path = "/home/thomas/Desktop/Master 2/Prog-Graphique-Appl-Indus/Prog-Graphique-et-Application-Industrielle./CSV/area.csv";
    QFile my_area(path);

    if(my_area.open(QFile::WriteOnly|QFile::Truncate)){
        QTextStream stream(&my_area);
        stream << "Aires," << "nb faces\n";
        for (auto& x: area_freq) {
            qDebug() << x.first << "," << x.second;
            if(x.second > 0)
                stream << x.first << "," << x.second << "\n";
        }
    }

    my_area.close();
    qDebug() << "CSV AREA WRITED !!";

    /** Angle dihedres **/
    //path =  "/Users/eliasmunoz/Documents/Git Projects/Prog-Graphique-et-Application-Industrielle./CSV/angledi.csv";
    //path =  "/home/kammerlocher/prog_indus/Prog-Graphique-et-Application-Industrielle./CSV/angleDihedre.csv";
    //path =  "/home/thomas/Desktop/Master 2/Prog-Graphique-Appl-Indus/Prog-Graphique-et-Application-Industrielle./CSV/angleDihedre.csv";
    path = "./CSV/angledi.csv";
    QFile my_dihedral(path);

    if(my_dihedral.open(QFile::WriteOnly|QFile::Truncate)){
        QTextStream stream(&my_dihedral);
        stream << "angles dihedres (degrés)," << "nb faces\n";
        for (auto& x: dihedral_freq) {
            qDebug() << x.first << "," << x.second;
            if(x.second > 0)
                stream << x.first << "," << x.second << "\n";
        }
    }

    /** Valences **/
    //path =  "/Users/eliasmunoz/Documents/Git Projects/Prog-Graphique-et-Application-Industrielle./CSV/angledi.csv";
    //path =  "/home/kammerlocher/prog_indus/Prog-Graphique-et-Application-Industrielle./CSV/valence.csv";
    path = "./CSV/valences.csv";
    QFile my_valence(path);

    if(my_valence.open(QFile::WriteOnly|QFile::Truncate)){
        QTextStream stream(&my_valence);
        stream << "Valences," << "Occurences\n";
        for (auto& x: valence_freq) {
            qDebug() << x.first << "," << x.second;
            if(x.second > 0)
                stream << x.first << "," << x.second << "\n";
        }
    }
    my_valence.close();
    /** Ecart Angulaire **/
    //path =  "/Users/eliasmunoz/Documents/Git Projects/Prog-Graphique-et-Application-Industrielle./CSV/angledi.csv";
    //path = "/home/thomas/Desktop/Master 2/Prog-Graphique-Appl-Indus/Prog-Graphique-et-Application-Industrielle./CSV/ecartAngulaire.csv";
    path = "./CSV/Ecart angulaire.csv";
    QFile my_ecart(path);

    if(my_ecart.open(QFile::WriteOnly|QFile::Truncate)){
        QTextStream stream(&my_ecart);
        stream << "Ecart angulaire," << "nb sommets\n";
        for(auto& x: ecart_angulaire) {
            qDebug() << x.first << "," << x.second;
            if(x.second > 0) stream << x.first << "," << x.second << "\n";
        }
    }
}

void createBox(MyMesh::Point min, MyMesh::Point max, MyMesh * _mesh){

    MyMesh::VertexHandle vhandle[8];
    qDebug() << "min : (" << min[0] << ", " << min[1] << ", " << min[2] << ")" << "\n";
    qDebug() << "max : (" << max[0] << ", " << max[1] << ", " << max[2] << ")" << "\n";

    vhandle[0] = _mesh->add_vertex(MyMesh::Point(min[0], min[1], max[2])); // -1 -1 1
    vhandle[1] = _mesh->add_vertex(MyMesh::Point(max[0], min[1], max[2])); // 1 -1 1
    vhandle[2] = _mesh->add_vertex(max); // 1 1 1
    vhandle[3] = _mesh->add_vertex(MyMesh::Point(min[0], max[1], max[2])); // -1 1 1
    vhandle[4] = _mesh->add_vertex(min); // -1 -1 -1
    vhandle[5] = _mesh->add_vertex(MyMesh::Point(max[0], min[1], min[2])); // 1 -1 -1
    vhandle[6] = _mesh->add_vertex(MyMesh::Point(max[0], max[1], min[2])); // 1 1 -1
    vhandle[7] = _mesh->add_vertex(MyMesh::Point(min[0], max[1], min[2])); // -1 1 -1

    std::vector<MyMesh::VertexHandle>  face_vhandles;

    // face 0 1 2
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[2]);
    _mesh->add_face(face_vhandles);

    // face 2 3 0
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[3]);
    face_vhandles.push_back(vhandle[0]);
    _mesh->add_face(face_vhandles);

    // face 7 6 5
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[5]);
    _mesh->add_face(face_vhandles);

    // face 5 4 7
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[4]);
    face_vhandles.push_back(vhandle[7]);
    _mesh->add_face(face_vhandles);

    // face 1 0 4
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[4]);
    _mesh->add_face(face_vhandles);

    // face 4 5 1
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[4]);
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[1]);
    _mesh->add_face(face_vhandles);

    // face 2 1 5
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[5]);
    _mesh->add_face(face_vhandles);

    // face 5 6 2
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[2]);
    _mesh->add_face(face_vhandles);

    // face 3 2 6
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[3]);
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[6]);
    _mesh->add_face(face_vhandles);

    // face 6 7 3
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[3]);
    _mesh->add_face(face_vhandles);

    // face 0 3 7
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[3]);
    face_vhandles.push_back(vhandle[7]);
    _mesh->add_face(face_vhandles);

    // face 7 4 0
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[4]);
    face_vhandles.push_back(vhandle[0]);
    _mesh->add_face(face_vhandles);
}

void MainWindow::boite_englobante(MyMesh* _mesh)
{
    MyMesh::Point  max_coord = MyMesh::Point(0, 0, 0);
    MyMesh::Point  min_coord = MyMesh::Point(100000000, 10000000, 10000000);
    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it)
    {
        VertexHandle vh = *v_it;
        MyMesh::Point p = _mesh->point(vh);
        if(max_coord[0] < p[0]){
            max_coord[0] = p[0];
        }
        if(max_coord[1] < p[1]){
            max_coord[1] = p[1];
        }
        if(max_coord[2] < p[2]){
            max_coord[2] = p[2];
        }
        if(min_coord[0] > p[0]){
            min_coord[0] = p[0];
        }
        if(min_coord[1] > p[1]){
            min_coord[1] = p[1];
        }
        if(min_coord[2] > p[2]){
            min_coord[2] = p[2];
        }
    }
    n_v_old = _mesh->n_vertices();
    n_f_old = _mesh->n_faces();
    n_e_old = _mesh->n_edges();

    qDebug() << "nb faces avant : " << n_f_old;
    createBox(min_coord, max_coord, _mesh);
    qDebug() << "nb faces apres : " << _mesh->n_faces();

    for(int i = n_e_old; i< _mesh->n_edges(); i++){
        EdgeHandle eh = _mesh->edge_handle(i);
        _mesh->set_color(eh, MyMesh::Color(255, 0, 0));
        _mesh->data(eh).thickness += 5;
    }
    for(int i = n_f_old; i< _mesh->n_faces(); i++){
        FaceHandle fh = _mesh->face_handle(i);
        _mesh->set_color(fh, MyMesh::Color(150, 150, 150));
    }
    /*for(MyMesh::EdgeIter curEdge = box.edges_begin(); curEdge != box.edges_end(); curEdge++){
        EdgeHandle eh = curEdge;
        box.set_color(eh, MyMesh::Color(255, 0, 0));
        box.data(eh).thickness += 5;
    }

    for(MyMesh::FaceIter curFace = box.faces_begin(); curFace != box.faces_end(); curFace++){
        FaceHandle fh = curFace;
        box.set_color(fh, MyMesh::Color(150, 150, 150));
    }*/
}

MyMesh::Point MainWindow::centre_gravite(MyMesh *_mesh){
    MyMesh::Point centre_grav(0.0,0.0,0.0);
    MyMesh::Point p;
    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it){
        p = _mesh->point(*v_it);
        centre_grav += p;
    }
    centre_grav /= _mesh->n_vertices();
    qDebug() << "Centre de gravité 1 : " << centre_grav[0] << " " << centre_grav[1] << " " << centre_grav[2];
    return centre_grav;
}

std::map<uint,int> MainWindow::valence(MyMesh* _mesh)
{
    int nb_sommets = _mesh->n_vertices();
    uint valences[nb_sommets];
    for (int i = 0; i<nb_sommets; i++)
    {
        valences[i] = 0;
    }
    int cpt = 0;
    //uint max = 0;

    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it)
    {
        VertexHandle vh = *v_it;
        valences[cpt] += _mesh->valence(vh);
        ++cpt;
    }

    std::map<uint, int> nb_sommets_valence; //nombre de sommets ayant la valence comme indice
    for(int i = 0; i < cpt ; i++)
    {
        nb_sommets_valence[valences[i]] = 0;
    }

    for(int i = 0; i < cpt ; i++)
    {
        nb_sommets_valence[valences[i]] += 1;
    }
    return  nb_sommets_valence;
}

double MainWindow::calcul_area(MyMesh::Point p[]){
    MyMesh::Point v1 = p[0];
    MyMesh::Point v2 = p[1];
    MyMesh::Point v3 = p[2];

    // on recupere le vecteur v1v2
    MyMesh::Point u = v2 - v1;
    // on recupere le vecteur v1v3
    MyMesh::Point v = v3 - v1;

    double ux = u[0], uy = u[1], uz = u[2];
    double vx = v[0], vy = v[1], vz = v[2];

    // on calcul les determinant du produit vectoriel
    volatile double i = (uy*vz - vy*uz);
    volatile double j = (ux*vz - vx*uz);
    volatile double k = (ux*vy - vx*uy);

    volatile double area2 = i - j + k;
    volatile double area = abs(area2)/2.0;

    return area;

}
bool AreSame(double a, double b)
{
    return fabs(a - b) < DBL_EPSILON;
}

std::map<double, int> MainWindow::area_frequency(MyMesh* _mesh) {
    std::map<double, int> freq;
    // On recupere le nombre de faces
    int n_faces = _mesh->n_faces();
    // On crée un tableau qui contiendra l'aire des faces
    double aires[n_faces];
    double minArea = 1000000.0;
    double maxArea = 0.0;
    double current_area = 0.0;
    int cpt = 0;
    //initialisation
    for(int i = 0; i<n_faces; i++){
        aires[i] = 0;
    }

    MyMesh::Point points[3];
    int cpt_points = 0;
    // On va maintenant parcourir les faces du maillage et calculer l'aire associé
    for(MyMesh::FaceIter current_face = _mesh->faces_begin(); current_face != _mesh->faces_end(); current_face++){
        // On recupere les point de la face triangulé
        cpt_points= 0;
        for(MyMesh::FaceVertexIter current_vertex = _mesh->fv_iter(*current_face); current_vertex.is_valid(); current_vertex++){
            VertexHandle vh = *current_vertex;
            points[cpt_points++] = _mesh->point(vh);
        }
        current_area = calcul_area(points);
        aires[cpt++] = current_area;
        if(minArea > current_area) minArea = current_area;
        if(maxArea < current_area) maxArea = current_area;
    }

    volatile double sigma = 0.1*(maxArea - minArea);
    double current = minArea;
    qDebug() << "aire min : " << minArea;
    qDebug() << "aire max : " << maxArea;

    if(minArea == maxArea){
        for(int i = 0; i< n_faces; i++){
            qDebug() << "aire = " << aires[i];
            freq[current] +=1;
        }
    }
    else{
        while(current <= maxArea){
            freq[current] = 0;
            for(int i = 0; i< n_faces; i++){
                if(is_in_range(aires[i], current, sigma)){
                    qDebug() << "freq[" << current << "] = " << freq[current];
                    freq[current] +=1;
                }
            }
            if(sigma == 0)
                break;
            current += sigma;
        }
    }
    return freq;
}

bool MainWindow::is_in_range(double valueTest, double a, double marginOfError){
    if(valueTest > a - marginOfError && valueTest <= a + marginOfError)
        return true;
    return false;
}

std::map<MyMesh::Scalar, int> MainWindow::dihedral_angles(MyMesh *_mesh){
    MyStats<MyMesh::Scalar> angles;
    std::map<MyMesh::Scalar, int> frequency;
    MyMesh::Scalar pi= M_PI;

    //Initialisation de la map
    int i = 0;
    while(i <= 360){
        frequency[i] = 0;
        i += 10;
    }

    // On recupere la valeur des angles diedre
    for(MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
        EdgeHandle eh = *curEdge;
        if(!_mesh->is_boundary(eh)){
            MyMesh::Scalar s_rad = _mesh->calc_dihedral_angle(eh);
            MyMesh::Scalar s_deg = (s_rad*180)/pi;
            angles.push_back(s_deg);
        }


    }
    // On On enumere le nombre d'angle pour chaque tranche de 10° de 0° a 360°
    std::vector<MyMesh::Scalar> v = angles.get_vector();
    for (int i = 0; i<=360 ; i+=10) {
        for (int j=0; j < v.size(); j++) {
            if(is_in_range(v.at(j), i, 5))
                frequency[i]++;
        }
    }

    /*for (int i = 0; i<=360 ; i+=10) {
        qDebug() << "Nombres d'angles pour" << i <<"(degres)"<<frequency[i];
    }*/

    return frequency;
}

//test si il y a la présence d'une face isolée.
bool MainWindow::test_lonely_face(MyMesh* _mesh){
    bool face_seule = true;

    for(MyMesh::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it){
        face_seule = true;
        for(MyMesh::FaceFaceIter ff_it = _mesh->ff_iter(*f_it); ff_it.is_valid(); ++ff_it) {
            face_seule = false;
        }
    }
    return face_seule;
}

//test si il y a la présence d'un points isolé
bool MainWindow::test_lonely_vertex(MyMesh* _mesh){
    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it){
        VertexHandle current_vertex = *v_it;
        if(_mesh->is_isolated(current_vertex)) return true;
    }
    return false;
}

//test si il y a la présence d'une arete isolé.
bool MainWindow::test_lonely_edge(MyMesh* _mesh){
    for(MyMesh::EdgeIter e_it = _mesh->edges_begin(); e_it != _mesh->edges_end(); ++e_it){
        EdgeHandle current_edge = *e_it;
        HalfedgeHandle he1 = _mesh->halfedge_handle(current_edge, 0);
        HalfedgeHandle he2 = _mesh->halfedge_handle(current_edge, 1);
        if(!_mesh->face_handle(he1).is_valid() && !_mesh->face_handle(he2).is_valid()) return true;
    }
    return false;
}

//Calcule la normale d'une face
MyMesh::Normal MainWindow::face_norm(MyMesh *_mesh, MyMesh::FaceHandle face){
    return _mesh->calc_face_normal(face);
}

// Affiche la normale de toute les faces
void MainWindow::print_faces_norm(MyMesh *_mesh){
    MyMesh::Normal normal;
    for(MyMesh::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it){
        normal = face_norm(_mesh, *f_it);
        qDebug() << "Normale face n°" << f_it->idx() << " : " << normal[0] << "," << normal[1] << "," << normal[2];
    }
}

//Renvoit la normale d'un point.
MyMesh::Normal MainWindow::vertex_norm(MyMesh *_mesh, MyMesh::VertexHandle vertex){
    return _mesh->calc_vertex_normal(vertex);
}

//affiche toutes les normales de tous les points.
void MainWindow::print_vertices_norm(MyMesh *_mesh){
    MyMesh::Normal normale;
    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it){
        normale = vertex_norm(_mesh, *v_it);
        qDebug() << "Normale point n°" << v_it->idx() << " : " << normale[0] << "," << normale[1] << "," << normale[2];
    }
}

//Calcule l'ecart angulaire pour chaque point.
void MainWindow::ecart_angulaire(MyMesh* _mesh){
    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it) {
        float current_angle = 0.0;
        MyMesh::Normal face_normal;
        MyMesh::Normal vertex_normal = _mesh->calc_vertex_normal(*v_it);
        for(MyMesh::VertexFaceIter vf_it = _mesh->vf_iter(*v_it); vf_it.is_valid(); ++vf_it){
            face_normal = _mesh->calc_face_normal(*vf_it);
            float norm_Vertex_Normal = sqrt(pow(vertex_normal[0],2) + pow(vertex_normal[1],2) + pow(vertex_normal[2],2));
            float norm_Face_Normal = sqrt(pow(face_normal[0],2) + pow(face_normal[1],2) + pow(face_normal[2],2));
            float prod_scalaire = dot(vertex_normal,face_normal);
            float new_angle = acos(prod_scalaire/(norm_Vertex_Normal*norm_Face_Normal));
            if(new_angle > current_angle) current_angle = new_angle;
        }
        _mesh->data(*v_it).value = current_angle;
    }
}

std::map<MyMesh::Scalar, int> MainWindow::ecart_ang(MyMesh* _mesh){
    qDebug() << __FUNCTION__;
    std::vector<MyMesh::Scalar> ecart_sommet;
    std::map<MyMesh::Scalar, int> ecart_ang;
    bool found_in_map = false;

    // On recupere la valeur de l'ecart angulaire pour chaque sommet.
    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it) {
        MyMesh::Scalar current_angle = 0.0;
        MyMesh::Normal face_normal;
        MyMesh::Normal vertex_normal = _mesh->calc_vertex_normal(*v_it);
        for(MyMesh::VertexFaceIter vf_it = _mesh->vf_iter(*v_it); vf_it.is_valid(); ++vf_it){
            face_normal = _mesh->calc_face_normal(*vf_it);
            float norm_Vertex_Normal = sqrt(pow(vertex_normal[0],2) + pow(vertex_normal[1],2) + pow(vertex_normal[2],2));
            float norm_Face_Normal = sqrt(pow(face_normal[0],2) + pow(face_normal[1],2) + pow(face_normal[2],2));
            float prod_scalaire = dot(vertex_normal,face_normal);
            float new_angle = acos(prod_scalaire/(norm_Vertex_Normal*norm_Face_Normal));
            if(new_angle > current_angle) current_angle = new_angle;
        }


//        _mesh->data(*v_it).value = current_angle;
        ecart_sommet.push_back(current_angle);
    }

    for(int i = 0; i < (int)_mesh->n_vertices(); i++){
        MyMesh::Scalar cle = (int)(ecart_sommet[i]*180/M_PI);
        ecart_ang[cle] = 0;
    }

    for(int i = 0; i < (int)_mesh->n_vertices(); i++){
        MyMesh::Scalar cle = (int)(ecart_sommet[i]*180/M_PI);
        ecart_ang[cle] += 1;
    }
    return ecart_ang;
}

int MainWindow::nb_faces_isole(MyMesh* _mesh){
    bool face_seule = true;
    int nb_faces_seules = 0;

    for(MyMesh::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it){
        face_seule = true;
        for(MyMesh::FaceFaceIter ff_it = _mesh->ff_iter(*f_it); ff_it.is_valid(); ++ff_it) {
            face_seule = false;
        }
        if(face_seule) nb_faces_seules++;
    }
    return nb_faces_seules;
}

int MainWindow::nb_points_isole(MyMesh* _mesh){
    int nb_points_seuls = 0;

    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it){
        VertexHandle current_vertex = *v_it;
        if(_mesh->is_isolated(current_vertex)) nb_points_seuls++;
    }

    return nb_points_seuls;
}

int MainWindow::nb_aretes_isole(MyMesh* _mesh){
    int nb_aretes_seules = 0;

    for(MyMesh::EdgeIter e_it = _mesh->edges_begin(); e_it != _mesh->edges_end(); ++e_it){
        EdgeHandle current_edge = *e_it;
        HalfedgeHandle he1 = _mesh->halfedge_handle(current_edge, 0);
        HalfedgeHandle he2 = _mesh->halfedge_handle(current_edge, 1);
        if(!_mesh->face_handle(he1).is_valid() && !_mesh->face_handle(he2).is_valid()) nb_aretes_seules++;
    }
    return nb_aretes_seules;
}

float MainWindow::Ai(MyMesh* _mesh, int vertexId){
    float aire = 0;

    VertexHandle currentvertex = _mesh->vertex_handle(vertexId);

    for(MyMesh::VertexFaceIter vf_it = _mesh->vf_iter(currentvertex); vf_it.is_valid(); ++vf_it){
        FaceHandle fh = *vf_it;
        HalfedgeHandle he = _mesh->halfedge_handle(fh);
        aire += _mesh->calc_sector_area(he);
    }

    return aire/3;
}

float MainWindow::angleFF(MyMesh *_mesh, int faceID0, int faceID1, int vertID0, int vertID1){
    float angle;

    _mesh->request_vertex_status();
    _mesh->request_edge_status();
    _mesh->request_face_status();

    /*get the faces*/
    MyMesh::FaceHandle *f0 = new FaceHandle(faceID0);
    MyMesh::FaceHandle *f1 = new FaceHandle(faceID1);

    /*calcul the normal faces*/
    MyMesh::Normal n0 = _mesh->calc_face_normal(*f0);
    MyMesh::Normal n1 = _mesh->calc_face_normal(*f1);
    float Norm_n0 = sqrt(pow(n0[0],2) + pow(n0[1],2) + pow(n0[2],2));
    float Norm_n1 = sqrt(pow(n1[0],2) + pow(n1[1],2) + pow(n1[2],2));
    float dot_n0_n1 = dot(n0,n1);

    /*u.v = ||u|| . ||v|| . cos(alpha)*/
    angle = acos(dot_n0_n1/(Norm_n0*Norm_n1));

    /*get the vertices*/
    VertexHandle *v0 = new VertexHandle(vertID0);
    VertexHandle *v1 = new VertexHandle(vertID1);

    /*make points from vertices to get coords (x, y, z)*/
    MyMesh::Point point0 = _mesh->point(*v0);
    MyMesh::Point point1 = _mesh->point(*v1);

    /* calcul the vector point0 -> point1 */
    double *vector =new double(3);
    vector[0] = point1[0] - point0[0];
    vector[1] = point1[1] - point0[1];
    vector[2] = point1[2] - point0[2];

    /*cross product between the two normals*/
    double *cross_product = new double(3);
    cross_product[0] = n0[1]*n1[2] - n0[2]*n1[1];
    cross_product[1] = n0[2]*n1[0] - n0[0]*n1[2];
    cross_product[2] = n0[0]*n1[1] - n0[1]*n1[0];

    float cpx = cross_product[0];
    float cpy = cross_product[1];
    float cpz = cross_product[2];

    float vx = vector[0];
    float vy = vector[1];
    float vz = vector[2];

    /*dot product entre les deux vecteurs pour avoir le signe de l'angle*/
    float dot_prod = cpx*vx + cpy*vy + cpz*vz;

    if(0 < dot_prod)
        return angle;
    else
        return -angle;
}


float MainWindow::angleEE(MyMesh* _mesh, int vertexID, int faceID){
    //Seulement besoin des deux points opposé a celui qu'on donne en entre de la fonction.
        float *x = new float[2];
        float *y = new float[2];
        float *z = new float[2];
        int cmpt=0;

        MyMesh::FaceHandle *f = new FaceHandle(faceID);

        VertexHandle baseVertex = _mesh->vertex_handle(vertexID);

        for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(*f); curVert.is_valid(); curVert ++)
        {

            VertexHandle vh = *curVert;


            if( vh == _mesh->vertex_handle(vertexID)) continue;
            else{
                x[cmpt] = _mesh->point(vh)[0];
                y[cmpt] = _mesh->point(vh)[1];
                z[cmpt] = _mesh->point(vh)[2];
                cmpt++;
            }
        }

        float *v1 = new float[3];
        float *v2 = new float[3];

        v1[0] = x[0] - _mesh->point(baseVertex)[0];
        v1[1] = y[0] - _mesh->point(baseVertex)[1];
        v1[2] = z[0] - _mesh->point(baseVertex)[2];

        v2[0] = x[1] - _mesh->point(baseVertex)[0];
        v2[1] = y[1] - _mesh->point(baseVertex)[1];
        v2[2] = z[1] - _mesh->point(baseVertex)[2];

        float prod = 0;
        for(int i = 0 ; i<3; i++){
            prod += v1[i] * v2[i];
        }
        float norv1 = sqrt(pow(v1[0],2)+pow(v1[1],2)+pow(v1[2],2));
        float norv2 = sqrt(pow(v2[0],2)+pow(v2[1],2)+pow(v2[2],2));

        float costet = prod/(norv1*norv2);
        float angle = acos(costet);

        return abs(angle);
}

void MainWindow::H_Curv(MyMesh* _mesh)
{
    float hCurv;
    float aire_bari;
    float somme;
    QList<float> longueurs;
    QList<float> angles;

    for(MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); v_it++){

        VertexHandle currentVertex = *v_it;
        hCurv=0;
        aire_bari=0;

        aire_bari = Ai(_mesh, currentVertex.idx());

        for(MyMesh::VertexVertexCWIter vv_cwiter = _mesh->vv_cwiter(currentVertex); vv_cwiter.is_valid(); vv_cwiter++){

            VertexHandle vertex = *vv_cwiter;
            HalfedgeHandle heh = _mesh->halfedge_handle(vertex);
            longueurs.push_back(_mesh->calc_edge_length(heh));

            FaceHandle f0 = _mesh->face_handle(heh);
            FaceHandle f1 = _mesh->face_handle(_mesh->opposite_halfedge_handle(heh));
            angles.push_back(angleFF(_mesh, f0.idx(), f1.idx(), currentVertex.idx(), vertex.idx()));

        }

        somme = 0;
        for(int i=0; i<longueurs.length(); i++){
            somme += angles[i] * longueurs[i];
        }

        float x = 1/(4*aire_bari);
        hCurv = x * somme;

        longueurs.clear();
        angles.clear();

        _mesh->data(currentVertex).value = hCurv;
    }
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    float kCurv;
    float angle;

    for (MyMesh::VertexIter v_it = _mesh->vertices_begin(); v_it!=_mesh->vertices_end(); v_it++)
        {
            VertexHandle vh = *v_it;
            kCurv = 0;
            angle = 0;
            for(MyMesh::VertexFaceIter  vf_it = _mesh->vf_iter(vh); vf_it; ++vf_it) {
                    FaceHandle fh = *vf_it;
                    angle += angleEE(_mesh, vh.idx(), fh.idx());
            }

            kCurv = (1/Ai(_mesh, vh.idx())) * ( 2*M_PI -  angle);
            _mesh->data(vh).value = kCurv;
        }
}

/* **** fin de la partie à compléter **** */

/* **** début de la partie boutons et IHM **** */
void MainWindow::on_pushButton_Affiche_Carac_clicked(){
//   qDebug() << __FUNCTION__;
   get_carac(&mesh);
}

void MainWindow::on_pushButton_normales_faces_clicked(){
//    qDebug() << __FUNCTION__;
//    print_faces_norm(&mesh);
    centre_gravite(&mesh);
}

void MainWindow::on_pushButton_normales_points_clicked(){
//    qDebug() << __FUNCTION__;
    print_vertices_norm(&mesh);
}

void MainWindow::on_pushButton_seuls_clicked(){
//    qDebug() << __FUNCTION__;
    qDebug() << "Face seule : " << test_lonely_face(&mesh) << ": il y a : " << nb_faces_isole(&mesh) << " faces seules";
    qDebug() << "Point seul : " << test_lonely_vertex(&mesh) << ": il y a : " << nb_points_isole(&mesh) << " points seuls";
    qDebug() << "Arete seule : " << test_lonely_edge(&mesh) << ": il y a : " << nb_aretes_isole(&mesh) << " aretes seules";
}

void MainWindow::on_pushButton_ecart_angulaire_clicked(){
    ecart_angulaire(&mesh);
    displayMesh(&mesh, DisplayMode::TemperatureMap);
}

void MainWindow::on_pushButton_H_clicked()
{
    H_Curv(&mesh);
    displayMesh(&mesh, DisplayMode::TemperatureMap); // permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_K_clicked()
{
    K_Curv(&mesh);
    displayMesh(&mesh, DisplayMode::TemperatureMap); // permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

/*  Cette fonction est à utiliser UNIQUEMENT avec le fichier testAngleArea.obj
    Elle est appelée par le bouton "Test angles/aires"

    Elle permet de vérifier les fonctions faceArea, angleFF et angleEE.
    Elle doit afficher :

    Aire de la face 0 : 2
    Aire de la face 1 : 2
    Angle entre les faces 0 et 1 : 1.5708
    Angle entre les faces 1 et 0 : -1.5708
    AngleEE au sommet 1 sur la face 0 : 0.785398 */
void MainWindow::on_pushButton_angleArea_clicked()
{

//    qDebug() << "Angle entre les faces 0 et 1 :" << angleFF(&mesh, 0, 1, 1, 2);
//    qDebug() << "Angle entre les faces 1 et 0 :" << angleFF(&mesh, 1, 0, 1, 2);

}

void MainWindow::on_pushButton_chargement_clicked()
{
    showBox = false;
    box_created = false;


    // fenêtre de sélection des fichiers
    fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);
    // on affiche le maillage
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, DisplayMode mode)
{

    /*if(showBox == true){
        MyMesh box = boite_englobante(_mesh);
        for(MyMesh::VertexIter curVert = box.vertices_begin(); curVert != box.vertices_end(); curVert++){
            VertexHandle vh = curVert;
            _mesh->add_vertex(box.point(vh));
        }

        for(MyMesh::EdgeIter curEdge = box.edges_begin(); curEdge != box.edges_end(); curEdge++){
            EdgeHandle eh = curEdge;

        }
    }*/

    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(mode == DisplayMode::TemperatureMap)
    {
        QVector<float> values;
        for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
            values.append(fabs(_mesh->data(*curVert).value));
        qSort(values);

        float range = values.at(values.size()*0.8);

        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){

                triCols[3*i+0] = 255;
                triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
            }
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::Normal)
    {
         MyMesh::ConstFaceVertexIter fvIt;
        if(showBox == true){
            for(int j = 0; j < n_f_old; j++){
                MyMesh::FaceHandle fIt = _mesh->face_handle(j);
                fvIt = _mesh->cfv_iter(fIt);
                triCols[3*i+0] = _mesh->color(fIt)[0]; triCols[3*i+1] = _mesh->color(fIt)[1]; triCols[3*i+2] = _mesh->color(fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                triCols[3*i+0] = _mesh->color(fIt)[0]; triCols[3*i+1] = _mesh->color(fIt)[1]; triCols[3*i+2] = _mesh->color(fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                triCols[3*i+0] = _mesh->color(fIt)[0]; triCols[3*i+1] = _mesh->color(fIt)[1]; triCols[3*i+2] = _mesh->color(fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++;
            }
        }
        else {
            MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
            MyMesh::ConstFaceVertexIter fvIt;
            for (; fIt!=fEnd; ++fIt)
            {
                fvIt = _mesh->cfv_iter(*fIt);
                triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;
                i++; ++fvIt;
                triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
                triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++;
            }
        }
    }

    if(mode == DisplayMode::ColorShading)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;


    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    export_csv();
}

void MainWindow::delete_box(MyMesh * _mesh){
    _mesh->request_face_status();
    _mesh->request_edge_status();
    _mesh->request_vertex_status();

    for(int i=n_f_old; i < _mesh->n_faces(); i++){
        FaceHandle fh= _mesh->face_handle(i);
        _mesh->delete_face(fh);
    }

    for(int i=n_e_old; i < _mesh->n_edges(); i++){
        EdgeHandle eh = _mesh->edge_handle(i);
        _mesh->delete_edge(eh);
    }

    for(int i=n_v_old; i < _mesh->n_vertices(); i++){
        VertexHandle vh = _mesh->vertex_handle(i);
        _mesh->delete_vertex(vh);
    }

}

void MainWindow::on_BoundingBox_clicked()
{
    if(showBox == false){
        showBox = true;
        if(box_created == false){
            boite_englobante(&mesh);
            box_created = true;
        }
    }
    else{
        showBox = false;
        box_created = false;
        OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());
        mesh.update_normals();
        resetAllColorsAndThickness(&mesh);
    }
    displayMesh(&mesh);
}

void MainWindow::on_spinBox_valueChanged(int arg1)
{
    seuil = arg1;
}

void MainWindow::on_pushButton_2_clicked()
{
    translate_to_origin(&mesh);
    displayMesh(&mesh);
}
