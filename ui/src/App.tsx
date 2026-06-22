import nci60Data from '@/data/nci60.json'
import golubData from '@/data/golub.json'
import { Hero } from '@/components/Hero'
import { DatasetCard } from '@/components/DatasetCard'
import { Playground } from '@/components/Playground'
import { InstallCTA } from '@/components/InstallCTA'
import type { DatasetOption } from '@/components/DatasetSelector'
import type { DemoData } from '@/types/demo'

const nci60 = nci60Data as unknown as DemoData
const golub = golubData as unknown as DemoData

const datasets: DatasetOption[] = [
  {
    key: 'nci60',
    label: 'NCI-60',
    description: 'Cancer cell lines across multiple tissue types',
    data: nci60,
  },
  {
    key: 'golub',
    label: 'Golub',
    description: 'Leukemia gene expression — AML vs B-ALL vs T-ALL',
    data: golub,
  },
]

export default function App() {
  const defaultData = nci60

  return (
    <div className="min-h-screen">
      <Hero />
      <DatasetCard data={defaultData} />
      <Playground datasets={datasets} defaultDataset="nci60" />
      <InstallCTA />
      <footer className="border-t border-slate-800 py-8 text-center text-xs text-slate-500">
        evaluomeR · GPL-3 ·{' '}
        <a
          href="https://github.com/neobernad/evaluomeR"
          className="text-slate-400 hover:text-white"
        >
          GitHub
        </a>
      </footer>
    </div>
  )
}
