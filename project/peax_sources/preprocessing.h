/*!
 * @author Marianna D'Addario, Dominik Kopczynski
 * @e-mail {Marianna.Daddario, Dominik.Kopczynski}@tu-dortmund.de
 * Copyright (c) 2013 Marianna D'Addario, Dominik Kopczynski
 *
 * LICENSE:
 * PEAX is a free software for academic use only.
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef PREPROCESSING_H
#define PREPROCESSING_H

#include "ims_core.h"
#include "fft.h"
#include <math.h>
#include <deque>

#define FFTFORWARD 1
#define FFTBACKWARD 0
#define DEF_MAX_EM_LOOPS 15

#define preprocessing(func_name) void func_name(IMSMeasurement* measurement, pmap* pipeline_parameters)

typedef void (*preprocessing_p)(IMSMeasurement*, pmap*);

preprocessing(de_tailing);
preprocessing(baseline_correction);
preprocessing(just_positives);
preprocessing(no_preprocessing);
preprocessing(de_noising);
preprocessing(crop);
preprocessing(mixed_smoothing);

#endif /* PREPROCESSING */