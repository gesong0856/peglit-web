<template>
  <div class="peglit-container">
    <header class="header">
      <h1 class="title">pegLIT</h1>
      <p class="subtitle">Automatically identify non-interfering nucleotide linkers between a pegRNA and 3' motif.</p>
    </header>

    <div class="input-table-wrapper">
      <table class="input-table">
        <thead>
          <tr>
            <th>Spacer</th>
            <th>Scaffold</th>
            <th>Template</th>
            <th>PBS</th>
            <th>Linker Pattern</th>
            <th>Motif</th>
          </tr>
        </thead>
        <tbody>
          <tr v-for="(row, idx) in inputRows" :key="idx">
            <td><input v-model="row.spacer" class="seq-input" /></td>
            <td><input v-model="row.scaffold" class="seq-input" /></td>
            <td><input v-model="row.template" class="seq-input" /></td>
            <td><input v-model="row.pbs" class="seq-input" /></td>

            <td class="linker-cell">
              <div v-if="loading" class="progress-container">
                <span class="countdown">{{ eta }}s</span>
                <div class="progress-bar">
                  <div class="progress-fill" :style="{ width: progress + '%' }"></div>
                </div>
              </div>
              <div v-else class="output-container">
                <input v-model="row.linkerPattern" class="seq-input" readonly placeholder="Output">
                <span class="char-count">{{ row.linkerPattern.length }} bp</span>
              </div>
            </td>

            <td><input v-model="row.motif" class="seq-input" /></td>
          </tr>
        </tbody>
      </table>

      <div class="table-buttons">
        <button class="btn-icon" @click="addRow">⊕</button>
        <button class="btn-icon" @click="importCSV">
          ↑<input type="file" ref="csvInput" accept=".csv" @change="handleCSV" style="display:none;">
        </button>
      </div>
    </div>

    <div class="start-btn-wrapper">
      <button class="btn-start" @click="startCalc" :disabled="loading">
        {{ loading ? 'RUNNING...' : 'START' }}
      </button>
    </div>

    <div class="result-wrapper" v-if="results.length">
      <h3>Optimized Linker Results</h3>
      <table class="result-table">
        <thead>
          <tr><th>#</th><th>Linker Sequence</th><th>Score</th></tr>
        </thead>
        <tbody>
          <tr v-for="(item, i) in results" :key="i">
            <td>{{ i+1 }}</td>
            <td class="seq-cell">{{ item.seq }}</td>
            <td>{{ item.score.join(' / ') }}</td>
          </tr>
        </tbody>
      </table>
      <button class="btn-fill" @click="fillTop">Fill Top Linker</button>
    </div>

    <div class="error-alert" v-if="errorMsg">{{ errorMsg }}</div>
  </div>
</template>

<script setup>
import { ref, onUnmounted } from 'vue'
import axios from 'axios'

// Use relative base URL so Vite proxy handles API routing in dev
axios.defaults.baseURL = ''
axios.defaults.headers.common['Content-Type'] = 'application/json'

const inputRows = ref([{
  spacer: '',
  scaffold: 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC',
  template: '',
  pbs: '',
  linkerPattern: '',
  motif: 'CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA'
}])

const loading = ref(false)
const progress = ref(0)
const eta = ref(15)
const results = ref([])
const errorMsg = ref('')
const csvInput = ref(null)

let taskId = ''
let interval = null

const addRow = () => inputRows.value.push({
  spacer: '', scaffold: '', template: '', pbs: '', linkerPattern: '', motif: ''
})

const importCSV = () => csvInput.value.click()
const handleCSV = (e) => {
  const file = e.target.files[0]
  const reader = new FileReader()
  reader.onload = (ev) => {
    const lines = ev.target.result.split('\n').filter(i => i.trim())
    inputRows.value = lines.slice(1).map(line => {
      const [s, sc, t, p, lp, m] = line.split(',')
      return { spacer: s||'', scaffold: sc||'', template: t||'', pbs: p||'', linkerPattern: '', motif: m||'' }
    })
  }
  reader.readAsText(file)
}

const startCalc = async () => {
  const valid = inputRows.value.filter(r => r.spacer.trim() && r.scaffold.trim() && r.template.trim() && r.pbs.trim())
  if (!valid.length) {
    errorMsg.value = '请填写 Spacer/Scaffold/Template/PBS'
    return
  }
  errorMsg.value = ''
  loading.value = true
  progress.value = 0
  eta.value = 15
  results.value = []

  try {
    // 强制打印请求日志，确认前端发请求
    console.log('[FRONTEND] Sending POST request to /api/peglit/batch')
    const { data } = await axios.post('/api/peglit/batch', { rows: valid })
    console.log('[FRONTEND] Received task_id:', data.task_id)
    taskId = data.task_id
    interval = setInterval(checkProgress, 500)
  } catch (e) {
    console.error('[FRONTEND] POST request failed:', e)
    errorMsg.value = `启动失败: ${e.message}`
    loading.value = false
  }
}

const checkProgress = async () => {
  try {
    const { data } = await axios.get(`/api/peglit/progress/${taskId}`)
    progress.value = data.progress
    
    // 真实倒计时
    if (data.progress > 0 && data.start_time) {
      const elapsed = (Date.now() - data.start_time * 1000) / 1000
      const total = elapsed / (data.progress / 100)
      eta.value = Math.max(0, Math.round(total - elapsed))
    }

    if (data.status === 'completed') {
      results.value = data.result.linker_list.map((seq, i) => ({
        seq, score: data.result.score_list[i]
      }))
      inputRows.value[0].linkerPattern = results.value[0].seq
      stop()
    }
    if (data.status === 'failed') {
      errorMsg.value = `计算失败: ${data.error}`
      stop()
    }
  } catch (e) {
    console.error('[FRONTEND] Progress check failed:', e)
    errorMsg.value = `进度获取失败: ${e.message}`
    stop()
  }
}

const stop = () => {
  clearInterval(interval)
  loading.value = false
}

const fillTop = () => {
  if (results.value.length) inputRows.value[0].linkerPattern = results.value[0].seq
}

onUnmounted(stop)
</script>

<style scoped>
/* 完全保留你原来的样式，一丝不动 */
.peglit-container {
  max-width: 1200px;
  margin: 0 auto;
  padding: 30px 20px;
  font-family: sans-serif;
}
.header { text-align: center; margin-bottom: 30px; }
.title { font-size: 70px; margin: 0; }
.subtitle { font-size: 22px; color: #666; }

.input-table-wrapper {
  background: white;
  border-radius: 8px;
  padding: 20px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}
.input-table { width: 100%; border-collapse: collapse; }
.input-table th {
  text-align: left;
  padding: 12px 8px;
  font-size: 16px;
  border-bottom: 1px solid #eee;
}
.input-table td {
  padding: 8px;
  position: relative;
}
.seq-input {
  width: 100%;
  padding: 12px 14px;
  border: 1px solid #ddd;
  border-radius: 4px;
  font-size: 14px;
  box-sizing: border-box;
  font-family: monospace;
}

.linker-cell { position: relative; height: 48px; display: flex; align-items: center; }
.progress-container {
  display: flex;
  align-items: center;
  gap: 10px;
  width: 100%;
}
.countdown {
  font-size: 14px;
  color: #666;
  font-style: italic;
  min-width: 35px;
}
.progress-bar {
  flex: 1;
  height: 6px;
  background: #e0e0e0;
  border-radius: 3px;
  overflow: hidden;
}
.progress-fill {
  height: 100%;
  background: #3f51b5;
  transition: width 0.3s ease;
}
.output-container { position: relative; width: 100%; }
.char-count {
  position: absolute;
  right: 14px;
  top: 50%;
  transform: translateY(-50%);
  font-size: 12px;
  color: #888;
  background: white;
  padding: 0 4px;
}

.table-buttons {
  display: flex;
  gap: 10px;
  margin-top: 10px;
}
.btn-icon {
  width: 36px; height: 36px;
  border: 1px solid #ddd;
  background: white;
  border-radius: 4px;
  cursor: pointer;
}

.start-btn-wrapper { text-align: center; margin: 30px 0; }
.btn-start {
  background: #3f51b5;
  color: white;
  border: none;
  padding: 14px 40px;
  border-radius: 6px;
  font-size: 18px;
  cursor: pointer;
}
.btn-start:disabled { background: #9fa8da; }

.result-wrapper {
  background: white;
  border-radius: 8px;
  padding: 20px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}
.result-table { width: 100%; border-collapse: collapse; }
.result-table th, .result-table td {
  padding: 8px;
  border: 1px solid #ddd;
  text-align: left;
}
.result-table th { background: #3f51b5; color: white; }
.seq-cell { font-family: monospace; }
.btn-fill {
  background: #4caf50;
  color: white;
  border: none;
  padding: 8px 16px;
  border-radius: 4px;
  cursor: pointer;
  margin-top: 10px;
}

.error-alert {
  color: red;
  background: #ffebee;
  padding: 10px;
  border-radius: 4px;
  text-align: center;
}
</style>
