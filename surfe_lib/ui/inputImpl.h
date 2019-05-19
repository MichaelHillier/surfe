#ifndef INPUTIMPL_H
#define INPUTIMPL_H

#include <modelling_parameters.h>
#include <ui_input.h>

#include <QDoubleValidator>
#include <QFileDialog>

class InputImpl : public QDialog, private Ui::SurfeInputDialog {
	Q_OBJECT;
public:
	InputImpl(QWidget *parent = 0);
	~InputImpl() {}
	static UI_Parameters GetDialogParameters();
	//UI_Parameters get_parameters();
public slots:
	void set_interface_data_file();
	void set_planar_data_file();
	void set_tangent_data_file();
	void set_inequality_data_file();
private:
	UI_Parameters get_parameters();
};

#endif
