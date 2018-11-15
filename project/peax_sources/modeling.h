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

#ifndef MODELING_H
#define MODELING_H

#include "ims_core.h"

// Segmentation
#define MINHEIGHT_SEGMENTATION 2

// EM Variables
#define CONVERGENCE_THRESHOLD 1e-2
#define REPLY 10


#define modeling(func_name) IMSPeakList* func_name(IMSMeasurement* measurement, IMSPeakList* peak_list, pmap* pipeline_parameters)

typedef IMSPeakList* (*modeling_p)(IMSMeasurement*, IMSPeakList*, pmap*);

modeling(pme_modeling);
modeling(no_modeling);

#endif /* MODELING_H */