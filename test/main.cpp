#include <modelling_input.h>
#include <continuous_property.h>
#include <modeling_methods.h>
#include <inputImpl.h>

#include <vector>


int main(int argc, char* argv[]) {

// 	QApplication app(argc, argv);
// 
// 	InputImpl *dialog = new InputImpl;
// 	dialog->show();
// 	app.exec();
// 	UI_Parameters parameters = dialog->get_parameters();

	UI_Parameters params = InputImpl::GetDialogParameters();

	return 0;
}
