#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "hmmer.h"

#if defined(__x86_64__)
#include <xmmintrin.h>
#endif

typedef struct {
  ESL_ALPHABET *abc;
  P7_HMM      **hmms;
  int           count;
} rosella_models;

typedef struct {
  ESL_DSQ  *dsq;
  int64_t   used, room;
  int64_t  *start, *length;
  uint64_t *id;
  int64_t   count, slots;
} rosella_targets;

void rosella_init(void) { p7_FLogsumInit(); }

void rosella_models_free(rosella_models *models) {
  if (models == NULL) return;
  for (int at = 0; at < models->count; at++) p7_hmm_Destroy(models->hmms[at]);
  free(models->hmms);
  esl_alphabet_Destroy(models->abc);
  free(models);
}

rosella_models *rosella_models_read(const char *path, char *error, size_t room) {
  P7_HMMFILE     *file = NULL;
  P7_HMM         *hmm = NULL;
  rosella_models *models = calloc(1, sizeof(*models));
  int             slots = 0, status;
  char            reason[eslERRBUFSIZE];

  if (models == NULL) return NULL;
  if (p7_hmmfile_OpenNoDB(path, NULL, &file, reason) != eslOK) {
    snprintf(error, room, "%s", reason);
    rosella_models_free(models);
    return NULL;
  }
  while ((status = p7_hmmfile_Read(file, &models->abc, &hmm)) == eslOK) {
    if (models->count == slots) {
      slots = slots ? 2 * slots : 256;
      P7_HMM **grown = realloc(models->hmms, sizeof(P7_HMM *) * slots);
      if (grown == NULL) {
        status = eslEMEM;
        p7_hmm_Destroy(hmm);
        break;
      }
      models->hmms = grown;
    }
    models->hmms[models->count++] = hmm;
  }
  if (status != eslEOF) {
    snprintf(error, room, "%s", file->errbuf);
    rosella_models_free(models);
    models = NULL;
  }
  p7_hmmfile_Close(file);
  return models;
}

int rosella_models_count(const rosella_models *models) { return models->count; }

int rosella_models_length(const rosella_models *models, int model) { return models->hmms[model]->M; }

rosella_targets *rosella_targets_new(void) { return calloc(1, sizeof(rosella_targets)); }

void rosella_targets_clear(rosella_targets *targets) {
  targets->used = 0;
  targets->count = 0;
}

void rosella_targets_free(rosella_targets *targets) {
  if (targets == NULL) return;
  free(targets->dsq);
  free(targets->start);
  free(targets->length);
  free(targets->id);
  free(targets);
}

int64_t rosella_targets_residues(const rosella_targets *targets) { return targets->used; }

// The residues go through the same map the FASTA reader uses, so a held protein scores as the
// one HMMER reads from a file.
static int digitise(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int64_t *residues, const char *protein,
                    int64_t length) {
  *residues = 0;
  dsq[0] = eslDSQ_SENTINEL;
  return esl_abc_dsqcat_noalloc(abc->inmap, dsq, residues, protein, length);
}

static void whole(ESL_SQ *sq, int64_t residues) {
  sq->n = sq->end = sq->W = sq->L = residues;
  sq->start = 1;
  sq->C = 0;
}

int rosella_targets_add(rosella_targets *targets, const rosella_models *models, uint64_t id,
                        const char *protein, int64_t length) {
  int64_t residues;

  if (targets->used + length + 2 > targets->room) {
    int64_t room = targets->room ? targets->room : 1 << 20;
    while (targets->used + length + 2 > room) room *= 2;
    ESL_DSQ *grown = realloc(targets->dsq, room);
    if (grown == NULL) return eslEMEM;
    targets->dsq = grown;
    targets->room = room;
  }
  if (targets->count == targets->slots) {
    int64_t slots = targets->slots ? 2 * targets->slots : 1 << 14;
    int64_t *start = realloc(targets->start, sizeof(int64_t) * slots);
    if (start == NULL) return eslEMEM;
    targets->start = start;
    int64_t *lengths = realloc(targets->length, sizeof(int64_t) * slots);
    if (lengths == NULL) return eslEMEM;
    targets->length = lengths;
    uint64_t *ids = realloc(targets->id, sizeof(uint64_t) * slots);
    if (ids == NULL) return eslEMEM;
    targets->id = ids;
    targets->slots = slots;
  }
  int status = digitise(models->abc, targets->dsq + targets->used, &residues, protein, length);
  if (status != eslOK) return status;
  targets->start[targets->count] = targets->used;
  targets->length[targets->count] = residues;
  targets->id[targets->count] = id;
  targets->count++;
  targets->used += residues + 2;
  return eslOK;
}

void rosella_text_free(char *text) { free(text); }

// One model over every held protein, set up as hmmsearch sets up each query, so each thread can
// take a model of its own. HMMER's x86 setup flushes denormals, so the thread gets its mode back.
char *rosella_search(const rosella_models *models, int model, const rosella_targets *targets,
                     double floor, size_t *size) {
  P7_HMM      *hmm = models->hmms[model];
  P7_BG       *bg = p7_bg_Create(models->abc);
  P7_PROFILE  *gm = p7_profile_Create(hmm->M, models->abc);
  P7_OPROFILE *om = p7_oprofile_Create(hmm->M, models->abc);
  P7_TOPHITS  *th = p7_tophits_Create();
  ESL_SQ      *sq = esl_sq_CreateDigital(models->abc);
  P7_PIPELINE *pli = NULL;
  char        *text = NULL;
  char         name[32];
  FILE        *sink = NULL;
  int          failed = 1;
#if defined(__x86_64__)
  unsigned int mode = _mm_getcsr();
#endif

  impl_Init();
  if (bg == NULL || gm == NULL || om == NULL || th == NULL || sq == NULL) goto DONE;
  p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL);
  p7_oprofile_Convert(gm, om);
  if ((pli = p7_pipeline_Create(NULL, om->M, 100, FALSE, p7_SEARCH_SEQS)) == NULL) goto DONE;
  pli->T = pli->domT = floor;
  pli->by_E = pli->dom_by_E = FALSE;
  pli->Z = pli->domZ = 1.0;
  pli->Z_setby = pli->domZ_setby = p7_ZSETBY_OPTION;
  pli->show_alignments = FALSE;
  if (p7_pli_NewModel(pli, om, bg) != eslOK) goto DONE;

  for (int64_t at = 0; at < targets->count; at++) {
    int64_t n = targets->length[at];
    if (esl_sq_GrowTo(sq, n) != eslOK) goto DONE;
    memcpy(sq->dsq, targets->dsq + targets->start[at], n + 2);
    snprintf(name, sizeof(name), "%llu", (unsigned long long) targets->id[at]);
    esl_sq_SetName(sq, name);
    whole(sq, n);
    p7_pli_NewSeq(pli, sq);
    p7_bg_SetLength(bg, n);
    p7_oprofile_ReconfigLength(om, n);
    if (p7_Pipeline(pli, om, bg, sq, NULL, th) != eslOK) goto DONE;
    esl_sq_Reuse(sq);
    p7_pipeline_Reuse(pli);
  }
  p7_tophits_SortBySortkey(th);
  p7_tophits_Threshold(th, pli);
  if ((sink = open_memstream(&text, size)) == NULL) goto DONE;
  failed = p7_tophits_TabularDomains(sink, hmm->name, hmm->acc, th, pli, FALSE) != eslOK;
  fclose(sink);

DONE:
  if (failed) {
    free(text);
    text = NULL;
  }
  esl_sq_Destroy(sq);
  p7_tophits_Destroy(th);
  p7_pipeline_Destroy(pli);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
#if defined(__x86_64__)
  _mm_setcsr(mode);
#endif
  return text;
}

// Set up as hmmalign sets up its one model, so the A2M rows are the ones it would write.
char *rosella_align(const char *profile, size_t profile_size, const char *const *proteins,
                    const int64_t *lengths, int count, size_t *size) {
  P7_HMMFILE   *file = NULL;
  P7_HMM       *hmm = NULL;
  ESL_ALPHABET *abc = NULL;
  ESL_MSA      *msa = NULL;
  ESL_SQ      **sq = calloc(count, sizeof(ESL_SQ *));
  P7_TRACE    **tr = calloc(count, sizeof(P7_TRACE *));
  char         *text = NULL;
  char          name[32];
  FILE         *sink = NULL;
  int           failed = 1;

  if (sq == NULL || tr == NULL) goto DONE;
  if (p7_hmmfile_OpenBuffer(profile, (int) profile_size, &file) != eslOK) goto DONE;
  if (p7_hmmfile_Read(file, &abc, &hmm) != eslOK) goto DONE;
  for (int at = 0; at < count; at++) {
    int64_t residues;
    if ((sq[at] = esl_sq_CreateDigital(abc)) == NULL) goto DONE;
    if ((tr[at] = p7_trace_CreateWithPP()) == NULL) goto DONE;
    if (esl_sq_GrowTo(sq[at], lengths[at]) != eslOK) goto DONE;
    if (digitise(abc, sq[at]->dsq, &residues, proteins[at], lengths[at]) != eslOK) goto DONE;
    whole(sq[at], residues);
    snprintf(name, sizeof(name), "%d", at);
    esl_sq_SetName(sq[at], name);
  }
  if (p7_tracealign_computeTraces(hmm, sq, 0, count, tr) != eslOK) goto DONE;
  if (p7_tracealign_Seqs(sq, tr, count, hmm->M, p7_ALL_CONSENSUS_COLS, hmm, &msa) != eslOK) goto DONE;
  if ((sink = open_memstream(&text, size)) == NULL) goto DONE;
  failed = esl_msafile_Write(sink, msa, eslMSAFILE_A2M) != eslOK;
  fclose(sink);

DONE:
  if (failed) {
    free(text);
    text = NULL;
  }
  for (int at = 0; sq != NULL && at < count; at++) esl_sq_Destroy(sq[at]);
  for (int at = 0; tr != NULL && at < count; at++) p7_trace_Destroy(tr[at]);
  free(sq);
  free(tr);
  esl_msa_Destroy(msa);
  p7_hmm_Destroy(hmm);
  if (file != NULL) p7_hmmfile_Close(file);
  esl_alphabet_Destroy(abc);
  return text;
}
