#include <inputImpl.h>
#include <QStyleFactory>

InputImpl::InputImpl(QWidget *parent /*= 0*/)
{
	setupUi(this);

	single_surface_button->setChecked(true);
	interface_frame->setVisible(false);
	planar_frame->setVisible(false);
	tangent_frame->setVisible(false);
	inequality_frame->setVisible(false);
	uncertainty_frame->setVisible(false);


	QDoubleValidator *doubleValues = new QDoubleValidator(this);
	doubleValues->setBottom(0.0);
	doubleValues->setDecimals(2);
	shape_parameter_lineedit->setValidator(doubleValues);
	smoothing_amount_lineedit->setValidator(doubleValues);
	interface_uncertainty_lineedit->setValidator(doubleValues);
	angular_uncertainty_lineedit->setValidator(doubleValues);

	this->layout()->setSizeConstraint(QLayout::SetFixedSize);

	adjustSize();

	connect(interface_browse_button, SIGNAL(clicked()), SLOT(set_interface_data_file()));
	connect(planar_browse_button, SIGNAL(clicked()), SLOT(set_planar_data_file()));
	connect(tangent_browse_button, SIGNAL(clicked()), SLOT(set_tangent_data_file()));
	connect(inequality_browse_button, SIGNAL(clicked()), SLOT(set_inequality_data_file()));
	connect(ok_button, SIGNAL(clicked()), SLOT(accept()));
}

UI_Parameters InputImpl::GetDialogParameters()
{

	char *argv[] = { "Surfe Input Dialog", nullptr };
	int argc = (int)(sizeof(argv)/sizeof(argv[0])) - 1;
	QApplication app(argc,argv);

	QPalette dark_palette = QPalette();
	dark_palette.setColor(QPalette::Window, QColor(53, 53, 53));
	dark_palette.setColor(QPalette::WindowText, Qt::gray);
	dark_palette.setColor(QPalette::Base, QColor(25, 25, 25));
	dark_palette.setColor(QPalette::AlternateBase, QColor(53, 53, 53));
	dark_palette.setColor(QPalette::ToolTipBase, Qt::white);
	dark_palette.setColor(QPalette::ToolTipText, Qt::white);
	dark_palette.setColor(QPalette::Text, Qt::black);
	dark_palette.setColor(QPalette::Button, QColor(53, 53, 53));
	dark_palette.setColor(QPalette::ButtonText, Qt::black);
	dark_palette.setColor(QPalette::BrightText, Qt::red);
	dark_palette.setColor(QPalette::Link, QColor(42, 130, 218));
	dark_palette.setColor(QPalette::Highlight, QColor(42, 130, 218));
	dark_palette.setColor(QPalette::HighlightedText, Qt::black);
	// give style to dialog
	app.setPalette(dark_palette);

	InputImpl *dialog = new InputImpl;
	dialog->show();
	app.exec();
	UI_Parameters parameters = dialog->get_parameters();
	return parameters;
}

UI_Parameters InputImpl::get_parameters()
{
	UI_Parameters params;
	if (single_surface_button->isChecked())
		params.model_type = Parameter_Types::ModelType::Single_surface;
	else if (lajuanie_button->isChecked())
		params.model_type = Parameter_Types::ModelType::Lajaunie_approach;
	else if (vector_field_button->isChecked())
		params.model_type = Parameter_Types::ModelType::Vector_field;
	else if (stratigraphic_button->isChecked())
		params.model_type = Parameter_Types::ModelType::Stratigraphic_horizons;
	else
		params.model_type = Parameter_Types::ModelType::Continuous_property;
	params.use_interface = interface_checkbox->isChecked();
	params.use_planar = planar_checkbox->isChecked();
	params.use_tangent = tangent_checkbox->isChecked();
	params.use_inequality = inequality_checkbox->isChecked();
	QString basis_text = basis_type_combobox->currentText();
	if (basis_text == "r3")
		params.basis_type = Parameter_Types::RBF::Cubic;
	else if (basis_text == "WendlandC2")
		params.basis_type = Parameter_Types::RBF::WendlandC2;
	else if (basis_text == "r")
		params.basis_type = Parameter_Types::RBF::R;
	else if (basis_text == "Gaussian")
		params.basis_type = Parameter_Types::RBF::Gaussian;
	else if (basis_text == "Multiquadratics")
		params.basis_type = Parameter_Types::RBF::MQ;
	else if (basis_text == "Thin Plate Spline")
		params.basis_type = Parameter_Types::RBF::TPS;
	else
		params.basis_type = Parameter_Types::RBF::IMQ; // Inverse Multiquadratics
	params.shape_parameter = shape_parameter_lineedit->text().toDouble();
	params.polynomial_order = polynomial_order_spinbox->text().toInt();
	params.model_global_anisotropy = use_global_plunge_model_checkbox->isChecked();
	params.use_restricted_range = restricted_range_checkbox->isChecked();
	params.use_regression_smoothing = use_regression_smoothing_checkbox->isChecked();
	params.smoothing_amount = smoothing_amount_lineedit->text().toDouble();
	params.interface_uncertainty = interface_uncertainty_lineedit->text().toDouble();
	params.angular_uncertainty = angular_uncertainty_lineedit->text().toDouble();
	std::string interface_file_str = interface_file_text->text().toStdString();
	std::string planar_file_str = planar_file_text->text().toStdString();
	std::string tangent_file_str = tangent_file_text->text().toStdString();
	std::string inequality_file_str = inequality_file_text->text().toStdString();
	// do some c string gymnastics 
	if (!interface_file_str.empty())
	{
		auto interface_file_length = interface_file_str.size();
		char *temp_interface_filename = new char[interface_file_length + 1];
		strncpy_s(
			temp_interface_filename,
			interface_file_length + 1, 
			interface_file_str.c_str(), 
			interface_file_length
		);
		temp_interface_filename[interface_file_length] = '\0';
		params.interface_file = temp_interface_filename;
	}
	if (!planar_file_str.empty())
	{
		auto planar_file_length = planar_file_str.size();
		char *temp_planar_filename = new char[planar_file_length + 1];
		strncpy_s(
			temp_planar_filename,
			planar_file_length + 1,
			planar_file_str.c_str(),
			planar_file_length
		);
		temp_planar_filename[planar_file_length] = '\0';
		params.planar_file = temp_planar_filename;
	}
	if (!tangent_file_str.empty())
	{
		auto tangent_file_length = tangent_file_str.size();
		char *temp_tangent_filename = new char[tangent_file_length + 1];
		strncpy_s(
			temp_tangent_filename,
			tangent_file_length + 1,
			tangent_file_str.c_str(),
			tangent_file_length
		);
		temp_tangent_filename[tangent_file_length] = '\0';
		params.tangent_file = temp_tangent_filename;
	}
	if (!inequality_file_str.empty())
	{
		auto inequality_file_length = inequality_file_str.size();
		char *temp_inequality_filename = new char[inequality_file_length + 1];
		strncpy_s(
			temp_inequality_filename,
			inequality_file_length + 1,
			inequality_file_str.c_str(),
			inequality_file_length
		);
		temp_inequality_filename[inequality_file_length] = '\0';
		params.inequality_file = temp_inequality_filename;
	}
	return params;
}

void InputImpl::set_interface_data_file()
{
	QString filename = QFileDialog::getOpenFileName(
		this, tr("Select Interface Data file"), "", tr("CSV Files (*.csv)"));
	QFileInfo file_info(filename);
	if (file_info.exists())
		interface_file_text->setText(filename);
}

void InputImpl::set_planar_data_file()
{
	QString filename = QFileDialog::getOpenFileName(
		this, tr("Select Planar Data file"), "", tr("CSV Files (*.csv)"));
	QFileInfo file_info(filename);
	if (file_info.exists())
		planar_file_text->setText(filename);
}

void InputImpl::set_tangent_data_file()
{
	QString filename = QFileDialog::getOpenFileName(
		this, tr("Select Tangent Data file"), "", tr("CSV Files (*.csv)"));
	QFileInfo file_info(filename);
	if (file_info.exists())
		tangent_file_text->setText(filename);
}

void InputImpl::set_inequality_data_file()
{
	QString filename = QFileDialog::getOpenFileName(
		this, tr("Select Inequality Data file"), "", tr("CSV Files (*.csv)"));
	QFileInfo file_info(filename);
	if (file_info.exists())
		inequality_file_text->setText(filename);
}

