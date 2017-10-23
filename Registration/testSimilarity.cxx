#include <iostream>
#include <string>
#include <chrono>
#include "stdio.h"
using namespace std;

int main(int argc, char *argv[])
{
	//char optimizer[] = "QuasiNewton";
	//char optimizer[] = "RegularGradientDescent";

	char optimizer[50];
	if (argc > 1){
		switch (atoi(argv[1]))
		{
		case 1:
			strcpy(optimizer, "QuasiNewton");
			break;
		case 2:
			strcpy(optimizer, "RegularGradientDescent");
			break;
		default:
			std::cerr << "Wrong optimizer code" << std::endl;
			std::cerr << "1: QuasiNewton"<< std::endl;
			std::cerr << "2: RegularGradientDescent" << std::endl;
			return EXIT_FAILURE;
		}
		
	}
	else{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " OptimizerCode" << endl;
		std::cerr << "1: QuasiNewton" << std::endl;
		std::cerr << "2: RegularGradientDescent" << std::endl;
		return EXIT_FAILURE;
	}

	double translationFactor[] = { 0.2};
	double rescaleFactor[] = { 2 };
	int iterationsNumber[] = { 500 };
	double learningRate[] = { 2 };
	double minimumStepLength[] = { 0.0001 };
	for (int k = 0; k < sizeof(iterationsNumber) / sizeof(iterationsNumber[0]); k++)
	{
		for (int r = 0; r < sizeof(rescaleFactor) / sizeof(rescaleFactor[0]); r++){
			for (int l = 0; l < sizeof(minimumStepLength) / sizeof(minimumStepLength[0]); l++){
				for (int j = 0; j < sizeof(learningRate) / sizeof(learningRate[0]); j++){
					for (int i = 0; i < sizeof(translationFactor) / sizeof(translationFactor[0]); i++){
						char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108";
						char* outputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/Similarity_0108second";

						char fixedImageFile[500];
						sprintf_s(fixedImageFile, "%s/first_base_line_ear.nii.gz", inputfolder);
						sprintf_s(fixedImageFile, "%s/second_base_line_ear.nii.gz", inputfolder);

						//sprintf_s(fixedImageFile, "%s/dicom_ear_male_leftEar.nii", inputfolder);
						
						char movingImageFile[300];
						sprintf_s(movingImageFile, "%s/atlasVolume_unlabeled", inputfolder);
												
						char registeredImageFile[300];
						sprintf_s(registeredImageFile, "%s/NormCorrelation_%s_Scale1_%g_%g_IterationNum%d_LearningRate%g_minimumStepLength%g_RegistrationOutput.nii", outputfolder, optimizer, translationFactor[i], rescaleFactor[r], iterationsNumber[k], learningRate[j], minimumStepLength[l]);

						char differenceBeforeFile[300];
						sprintf_s(differenceBeforeFile, "%s/NormCorrelation_%s_Scale1_%g_%g_IterationNum%d_LearningRate%g_minimumStepLengthg%g_differenceBeforeRegistration.nii", outputfolder, optimizer, translationFactor[i], rescaleFactor[r], iterationsNumber[k], learningRate[j], minimumStepLength[l]);

						char differenceAfterFile[300];
						sprintf_s(differenceAfterFile, "%s/NormCorrelation_%s_Scale1_%g_%g_IterationNum%d_LearningRate%g_minimumStepLength%g_differenceAfterRegistration.nii", outputfolder, optimizer, translationFactor[i], rescaleFactor[r], iterationsNumber[k], learningRate[j], minimumStepLength[l]);

						char labelmapOutput_SL[350];
						sprintf_s(labelmapOutput_SL, " %s/NormCorrelation_%s_Scale1_%g_%g_IterationNum%d_LearningRate%g_minimumStepLength%g_SL_labelmap_registrationOutput.nii", outputfolder, optimizer, translationFactor[i], rescaleFactor[r], iterationsNumber[k], learningRate[j], minimumStepLength[l]);

						char labelmapOutput_SM[350];
						sprintf_s(labelmapOutput_SM, " %s/NormCorrelation_%s_Scale1_%g_%g_IterationNum%d_LearningRate%g_minimumStepLength%g_SM_labelmap_registrationOutput.nii", outputfolder, optimizer, translationFactor[i], rescaleFactor[r], iterationsNumber[k], learningRate[j], minimumStepLength[l]);

						char labelmapOutput_ST[350];
						sprintf_s(labelmapOutput_ST, " %s/NormCorrelation_%s_Scale1_%g_%g_IterationNum%d_LearningRate%g_minimumStepLength%g_ST_labelmap_registrationOutput.nii", outputfolder, optimizer, translationFactor[i], rescaleFactor[r], iterationsNumber[k], learningRate[j], minimumStepLength[l]);

						char labelmapOutput_SV[350];
						sprintf_s(labelmapOutput_SV, " %s/NormCorrelation_%s_Scale1_%g_%g_IterationNum%d_LearningRate%g_minimumStepLength%g_SV_labelmap_registrationOutput.nii", outputfolder, optimizer, translationFactor[i], rescaleFactor[r], iterationsNumber[k], learningRate[j], minimumStepLength[l]);


						char labelmapInput_SV[200];
						sprintf_s(labelmapInput_SV, " %s/atlasVolume_SV_labelmap.nii", inputfolder);
						char labelmapInput_ST[200];
						sprintf_s(labelmapInput_ST, " %s/atlasVolume_ST_labelmap.nii", inputfolder);
						char labelmapInput_SM[200];
						sprintf_s(labelmapInput_SM, " %s/atlasVolume_SM_labelmap.nii", inputfolder);
						char labelmapInput_SL[200];
						sprintf_s(labelmapInput_SL, " %s/atlasVolume_SL_labelmap.nii", inputfolder);
						char parameter[50];
						sprintf_s(parameter, " %g %g %d %g %g", translationFactor[i], rescaleFactor[r], iterationsNumber[k], learningRate[j], minimumStepLength[l]);

						char command[3000] = "";
						if (strcmp(optimizer,"QuasiNewton")==0){
							strcat_s(command, "D:/ITK4.8/ITK_Projects_build/Registration/Debug/RegistrationSimilarity3DQN.exe ");
						}
						else if (strcmp(optimizer, "RegularGradientDescent") == 0)
						{
							strcat_s(command, "D:/ITK4.8/ITK_Projects_build/Registration/Debug/RegistrationSimilarity3D.exe ");
						}
						else
						{
							cout << "Wrong optimizer name!" <<endl;
							return EXIT_SUCCESS;
						}

						strcat_s(command, fixedImageFile);
						strcat_s(command, movingImageFile);
						strcat_s(command, registeredImageFile);

						strcat_s(command, differenceBeforeFile);
						strcat_s(command, differenceAfterFile);
						strcat_s(command, parameter);

						strcat_s(command, labelmapInput_SV);
						strcat_s(command, labelmapOutput_SV);
						strcat_s(command, labelmapInput_ST);
						strcat_s(command, labelmapOutput_ST);
						strcat_s(command, labelmapInput_SM);
						strcat_s(command, labelmapOutput_SM);
						strcat_s(command, labelmapInput_SL);
						strcat_s(command, labelmapOutput_SL);
						cout << command << endl;
						std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
						
						system(command);
						std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

						std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes"<<std::endl;

					}
				}
			}
		}
	}
	//system("pause");
	return EXIT_SUCCESS;
	//char *rigisredImageFile = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/Rigid_NormCorrelation_RegularGradientDescent_changedAltlas_scale1_5_100/CenteredRigidRegistration.nii";
	//argv[4] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/Rigid_NormCorrelation_RegularGradientDescent_changedAltlas_scale1_5_100/differenceBeforeRegistration.nii";
	//argv[5] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/Rigid_NormCorrelation_RegularGradientDescent_changedAltlas_scale1_5_100/differenceAfterRegistration.nii";


	////argv[1] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/ITK/ImageVolume/dicom_ear_male_leftEar.nii";
	////argv[2] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/ITK/ImageVolume/atlasVolume_unlabeled.nii";
	////argv[3] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/ITK/ImageVolume/Rigid_NormCorrelation_RegularGradientDescent_changedAltlas_scale1_5_200/CenteredRigidRegistration.nii";
	////argv[4] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/ITK/ImageVolume/Rigid_NormCorrelation_RegularGradientDescent_changedAltlas_scale1_5_200/differenceBeforeRegistration.nii";
	////argv[5] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/ITK/ImageVolume/Rigid_NormCorrelation_RegularGradientDescent_changedAltlas_scale1_5_200/differenceAfterRegistration.nii";


	//argv[0] = "RegistrationRigid3D";
	//	argv[1] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Image Volume/dicom_ear_male_leftEar.nii";
	//	argv[2] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Image Volume/atlasVolume_unlabeled.nii";
	//	argv[3] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Image Volume/Similarity_1_5_10_500/RegistrationResult.nii";
	//	argv[4] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Image Volume/Similarity_1_5_10_500/differenceBeforeRegistration.nii";
	//	argv[5] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Image Volume/Similarity_1_5_10_500/differenceAfterRegistration.nii";


	//argv[1] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/Borkholder - Micropump/Image Volume/dicom_ear_male_leftEar.nii";
	//argv[2] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/Borkholder - Micropump/Image Volume/atlasVolume_unlabeled.nii";
	//argv[3] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/Borkholder - Micropump/Image Volume/Similarity_1_5_10_500/RegistrationResult.nii";
	//argv[4] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/Borkholder - Micropump/Image Volume/Similarity_1_5_10_500/differenceBeforeRegistration.nii";
	//argv[5] = "/Users/zhenlinxu/baiduyun/Dropbox/Research/Borkholder - Micropump/Image Volume/Similarity_1_5_10_500/differenceAfterRegistration.nii";
}