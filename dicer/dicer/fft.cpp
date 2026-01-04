#include "polynomial.h"

#include "mkl_dfti.h"

namespace spaceless::detail {

vec_complex fft(const vec_complex &input, std::span<const int64_t> sizes) {
	auto output = input;
	DFTI_DESCRIPTOR_HANDLE desc;
	MKL_LONG status;
	if (sizes.size() == 1)
		status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, sizes.front());
	else
		status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, sizes.size(), sizes.data());
	status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status = DftiCommitDescriptor(desc);
	status = DftiComputeForward(desc, const_cast<std::complex<double> *>(input.data()), output.data());
	status = DftiFreeDescriptor(&desc);
	return output;
}

vec_complex ifft(const vec_complex &input, std::span<const int64_t> sizes) {
	auto output = input;
	DFTI_DESCRIPTOR_HANDLE desc;
	MKL_LONG status;
	if (sizes.size() == 1)
		status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, sizes.front());
	else
		status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, sizes.size(), sizes.data());
	status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / input.size());
	status = DftiCommitDescriptor(desc);
	status = DftiComputeBackward(desc, const_cast<std::complex<double> *>(input.data()), output.data());
	status = DftiFreeDescriptor(&desc);
	return output;
}

vec_complex dot(const vec_complex &lhs, const vec_complex &rhs) {
	vec_complex ret = lhs;
	for (int i = 0; i < ret.size(); i++)
		ret[i] *= rhs[i];
	return ret;
}

}