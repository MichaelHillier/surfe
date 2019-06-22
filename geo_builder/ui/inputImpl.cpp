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

InputParameters InputImpl::GetDialogParameters()
{
	char *argv[] = { "Surfe Input Dialog", nullptr };
	int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
	QApplication app(argc, argv);

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
	InputParameters parameters = dialog->get_parameters();
	return parameters;
}

InputParameters InputImpl::get_parameters()
{
	InputParameters input;
	Parameters params;
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
	input.interface_file = interface_file_text->text().toStdString();
	input.planar_file = planar_file_text->text().toStdString();
	input.tangent_file = tangent_file_text->text().toStdString();
	input.inequality_file = inequality_file_text->text().toStdString();

	input.parameters = params;

	return input;
}

void InputImpl::set_interface_data_file()
{
	QString filename = QFileDialog::getOpenFileName(
		this, tr("Select Interface Data file"), "", tr("Surfe Files (*.csv *.vtk *.vtp)"));
	QFileInfo file_info(filename);
	if (file_info.exists())
		interface_file_text->setText(filename);
}

void InputImpl::set_planar_data_file()
{
	QString filename = QFileDialog::getOpenFileName(
		this, tr("Select Planar Data file"), "", tr("Surfe Files (*.csv *.vtk *.vtp)"));
	QFileInfo file_info(filename);
	if (file_info.exists())
		planar_file_text->setText(filename);
}

void InputImpl::set_tangent_data_file()
{
	QString filename = QFileDialog::getOpenFileName(
		this, tr("Select Tangent Data file"), "", tr("Surfe Files (*.csv *.vtk *.vtp)"));
	QFileInfo file_info(filename);
	if (file_info.exists())
		tangent_file_text->setText(filename);
}

void InputImpl::set_inequality_data_file()
{
	QString filename = QFileDialog::getOpenFileName(
		this, tr("Select Inequality Data file"), "", tr("Surfe Files (*.csv *.vtk *.vtp)"));
	QFileInfo file_info(filename);
	if (file_info.exists())
		inequality_file_text->setText(filename);
}