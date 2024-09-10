#include "polynomial.h"

#ifdef USE_ONEMKL
#include "mkl_dfti.h"
#else
#define POCKETFFT_NO_MULTITHREADING
#include "external/pocketfft/pocketfft_hdronly.h"
#undef POCKETFFT_NO_MULTITHREADING
#endif

namespace spaceless::detail {

vec_complex fft(const vec_complex &input, std::span<const int64_t> sizes) {
	auto output = input;
#ifdef USE_ONEMKL
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
#else
	pocketfft::shape_t shape(sizes, sizes + dim), axes(dim);
	pocketfft::stride_t stride(dim);
	for (int64_t i = shape.size() - 1, j = sizeof(std::complex<double>); i >= 0; i--) {
		stride[i] = j;
		j *= shape[i];
		axes[i] = i;
	}
	pocketfft::c2c(shape, stride, stride, axes, pocketfft::FORWARD, input.data(), output.data(), 1.0);
#endif
	return output;
}

vec_complex ifft(const vec_complex &input, std::span<const int64_t> sizes) {
	auto output = input;
#ifdef USE_ONEMKL
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
#else
	pocketfft::shape_t shape(sizes, sizes + dim), axes(dim);
	pocketfft::stride_t stride(dim);
	for (int64_t i = shape.size() - 1, j = sizeof(std::complex<double>); i >= 0; i--) {
		stride[i] = j;
		j *= shape[i];
		axes[i] = i;
	}
	pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD, input.data(), output.data(), 1.0 / input.size());
#endif
	return output;
}

vec_complex dot(const vec_complex &lhs, const vec_complex &rhs) {
	vec_complex ret = lhs;
	for (int i = 0; i < ret.size(); i++)
		ret[i] *= rhs[i];
	return ret;
}

}