#include <iostream>
#include <string>
#include "stdio.h"
#include <chrono>

using namespace std;

int main()
{
	double translationFactor[] = { 0.2,0.3 };
	int iterationsNumber[] = {500 };
	double learningRate[] = { 1 };
	double minimumStepLength[] = { 0.0001 };
	for (int j = 0; j < sizeof(learningRate) / sizeof(learningRate[0]); j++)
	{
		for (int k = 0; k < sizeof(iterationsNumber) / sizeof(iterationsNumber[0]); k++)
		{
			for (int i = 0; i < sizeof(translationFactor) / sizeof(translationFactor[0]); i++)
			{
				/*char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20150801";
				char* outputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20150801/Affine/Scale1_%g_IterationNum%d_LearningRate%g_";
				char fixedImageFile[300];
				sprintf_s(fixedImageFile, "%s/dicom_ear_male_leftEar.nii", inputfolder);*/

				char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";

				char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/Similarity_0108first/NormCorrelation_RegularGradientDescent_Scale1_0.3_2_IterationNum500_LearningRate1_minimumStepLength0.0001_";
				char outputfolder[300];
				sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/SimilarityAffine_0108second/Scale1_%g_IterationNum%d_LearningRate%g_", translationFactor[i], iterationsNumber[k], learningRate[j]);
				char movingImageFile[300];
				sprintf_s(movingImageFile, "%sRegistrationOutput.nii", inputfolder);

				//char fixedImageFile[300];
				//sprintf_s(fixedImageFile, "%sfirst_base_line_ear.nii.gz", fixedfolder);

				//char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";
				/*char outputfolder[300];
				sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/Affine_0108second/Scale1_%g_IterationNum%d_LearningRate%g_", translationFactor[i], iterationsNumber[k], learningRate[j]);*/
				char fixedImageFile[300];
				sprintf_s(fixedImageFile, "%ssecond_base_line_ear.nii.gz", fixedfolder);
				/*char movingImageFile[300];
				sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);*/
				

				char registeredImageFile[200];
				sprintf_s(registeredImageFile, " %sRegistrationOutput.nii", outputfolder);

				char differenceBeforeFile[200];
				sprintf_s(differenceBeforeFile, " %sdifferenceBeforeRegistration.nii", outputfolder );

				char differenceAfterFile[200];
				sprintf_s(differenceAfterFile, " %sdifferenceAfterRegistration.nii", outputfolder );

				char labelmapOutput_SL[200];
				sprintf_s(labelmapOutput_SL, " %sSL_labelmap_registrationOutput.nii", outputfolder );

				char labelmapOutput_SM[200];
				sprintf_s(labelmapOutput_SM, " %sSM_labelmap_registrationOutput.nii", outputfolder);

				char labelmapOutput_ST[200];
				sprintf_s(labelmapOutput_ST, " %sST_labelmap_registrationOutput.nii", outputfolder);

				char labelmapOutput_SV[200];
				sprintf_s(labelmapOutput_SV, " %sSV_labelmap_registrationOutput.nii", outputfolder);

				//char labelmapInput_SV[200];
				//sprintf_s(labelmapInput_SV, " %s/atlasVolume_SV_labelmap.nii", inputfolder);
				//char labelmapInput_ST[200];
				//sprintf_s(labelmapInput_ST, " %s/atlasVolume_ST_labelmap.nii", inputfolder);
				//char labelmapInput_SM[200];
				//sprintf_s(labelmapInput_SM, " %s/atlasVolume_SM_labelmap.nii", inputfolder);
				//char labelmapInput_SL[200];
				//sprintf_s(labelmapInput_SL, " %s/atlasVolume_SL_labelmap.nii", inputfolder);

				char labelmapInput_SV[300];
				sprintf_s(labelmapInput_SV, " %sSV_labelmap_registrationOutput.nii", inputfolder);
				char labelmapInput_ST[300];
				sprintf_s(labelmapInput_ST, " %sST_labelmap_registrationOutput.nii", inputfolder);
				char labelmapInput_SM[300];
				sprintf_s(labelmapInput_SM, " %sSM_labelmap_registrationOutput.nii", inputfolder);
				char labelmapInput_SL[300];
				sprintf_s(labelmapInput_SL, " %sSL_labelmap_registrationOutput.nii", inputfolder);

				char parameter[50];
				sprintf_s(parameter, " %g %d %g", learningRate[j], iterationsNumber[k], translationFactor[i] );

				char command[3000] = "D:/ITK4.8/ITK_Projects_build/Registration/Debug/RegistrationAffine3D.exe ";

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

				std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes" << std::endl;
			}
		}
	}
	return 0;
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