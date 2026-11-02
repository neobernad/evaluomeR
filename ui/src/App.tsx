import nci60Data from '@/data/nci60.json'
import golubData from '@/data/golub.json'
import type { DemoData } from '@/types/demo'
import type { Dataset } from '@/lib/demo'
import { Intro } from '@/components/Intro'
import { Explorer } from '@/components/Explorer'
import { Install } from '@/components/Install'

const golub = golubData as unknown as DemoData
const nci60 = nci60Data as unknown as DemoData

const datasets: Dataset[] = [
  { key: 'golub', label: 'Golub leukemia', short: 'Golub', noun: 'leukemia samples', data: golub },
  { key: 'nci60', label: 'NCI-60 cell lines', short: 'NCI-60', noun: 'cancer cell lines', data: nci60 },
]

export default function App() {
  return (
    <>
      <Intro data={golub} />
      <Explorer datasets={datasets} />
      <Install />
    </>
  )
}
